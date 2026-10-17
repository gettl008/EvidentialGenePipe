"""Stage 05 - call each transcript reference or novel.

A transcript is compared to the reference *in CDS space only*: the transcript's
ORF interval is projected through the genome alignment, and only those genomic
blocks are tested against reference CDS exons.

Classes
  reference          sense CDS-exon overlap with a reference coding gene
  novel_antisense    CDS blocks overlap a reference CDS exon on the other strand
  novel_intronic     inside a reference coding gene's span, no CDS-exon overlap
  novel_intergenic   mapped confidently to a region with no reference coding gene
  novel_unmapped     no alignment, or below the coverage / identity thresholds
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import pysam

from .. import refannot
from ..utils import LOG, require_tools, run

STAGE = "05_classify"

REF_CLASS = "reference"
NOVEL = ("novel_antisense", "novel_intronic", "novel_intergenic", "novel_unmapped")


# --------------------------------------------------------------------- blocks
def aligned_segments(rec) -> list[tuple[int, int, int, int]]:
    """Return [(qstart, qend, rstart, rend)] 0-based half-open, in stored-read orientation.

    One segment per run of match ops; N (and long D) are not merged across.
    """
    segs = []
    q = 0
    r = rec.reference_start
    cur = None
    for op, ln in rec.cigartuples or []:
        if op in (0, 7, 8):            # M, =, X
            if cur is None:
                cur = [q, q + ln, r, r + ln]
            else:
                cur[1] += ln
                cur[3] += ln
            q += ln
            r += ln
        elif op == 1:                  # I
            if cur is not None:
                cur[1] += ln
            q += ln
        elif op == 2:                  # D
            if cur is not None:
                cur[3] += ln
            r += ln
        elif op == 3:                  # N
            if cur is not None:
                segs.append(tuple(cur))
                cur = None
            r += ln
        elif op in (4,):               # S
            q += ln
        elif op in (5, 6):             # H, P
            pass
    if cur is not None:
        segs.append(tuple(cur))
    return segs


def aln_metrics(rec) -> tuple[int, int, float]:
    """(aligned_query_bases, query_length, identity)."""
    m = ins = dele = 0
    for op, ln in rec.cigartuples or []:
        if op in (0, 7, 8):
            m += ln
        elif op == 1:
            ins += ln
        elif op == 2:
            dele += ln
    qlen = rec.infer_read_length() or 0
    nm = rec.get_tag("NM") if rec.has_tag("NM") else 0
    mism = max(0, nm - ins - dele)
    ident = (m - mism) / m if m else 0.0
    return m + ins, qlen, ident


def project_interval(rec, qs: int, qe: int) -> list[tuple[int, int]]:
    """Project a stored-orientation query interval [qs,qe) onto genomic blocks."""
    out = []
    for a, b, ra, rb in aligned_segments(rec):
        # segment maps query [a,b) to reference [ra,rb) with internal indels;
        # approximate proportionally at the segment level (indel-scale error only)
        lo, hi = max(a, qs), min(b, qe)
        if hi <= lo:
            continue
        span_q = b - a
        span_r = rb - ra
        if span_q <= 0:
            continue
        s = ra + int(round((lo - a) * span_r / span_q))
        e = ra + int(round((hi - a) * span_r / span_q))
        if e > s:
            out.append((s, e))
    return out


def overlap_bp(a0, a1, b0, b1) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


# ------------------------------------------------------------------ ref links
def load_refcds_map(cfg, ra: refannot.RefAnnot) -> dict[str, str]:
    """refcds_NNNNNN -> reference gene ID."""
    idmap = cfg.outdir / "01_inputs" / "ref_cds.idmap.tsv"
    out: dict[str, str] = {}
    if not idmap.exists():
        return out
    with open(idmap) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = ra.tx2gene.get(f[1], f[1])
    LOG.info("reference CDS track: %d sequences linked to genes", len(out))
    return out


def load_dupcounts(cfg) -> dict[str, int]:
    p = cfg.outdir / "01_inputs" / "flnc.idmap.tsv"
    out: dict[str, int] = {}
    if not p.exists():
        return out
    with open(p) as fh:
        header = next(fh, "").rstrip("\n").split("\t")
        if "dupcount" not in header:
            return out
        i = header.index("dupcount")
        j = header.index("new_id")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > i:
                out[f[j]] = int(f[i])
    return out


def quantify_flnc(cfg, transcripts_fa: Path) -> dict[str, int]:
    """Map FLNC reads back to the final transcriptome; count primary alignments."""
    require_tools(["minimap2"])
    flnc = cfg.outdir / "01_inputs" / "flnc.fa"
    if not flnc.exists():
        return {}
    out = cfg.dir("05_classify")
    paf = out / f"{cfg.sample}.flnc_vs_transcripts.paf"
    if paf.exists() and paf.stat().st_size > 0 and not getattr(cfg, "force", False):
        LOG.info("reusing existing %s", paf.name)
    else:
        run(["minimap2", "-x", "map-hifi", "-t", str(cfg.threads), "--secondary=no",
             str(transcripts_fa), str(flnc)], stdout=paf)
    dup = load_dupcounts(cfg)
    counts: dict[str, int] = defaultdict(int)
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 12:
                continue
            # FLNC reads and evigene transcripts are both ORF-sense oriented, so a
            # reverse-strand hit is an antisense read, not support for this transcript
            if f[4] != "+":
                continue
            counts[f[5]] += dup.get(f[0], 1)
    LOG.info("FLNC support counted for %d transcripts", len(counts))
    return counts


# --------------------------------------------------------------------- driver
def run_stage(cfg, harvest: dict, bam_path: Path) -> Path:
    ra = refannot.load(cfg.annotation, cfg.d["inputs"].get("annotation_format", "auto"))
    refcds2gene = load_refcds_map(cfg, ra)
    cl = cfg.d["classify"]
    min_cov = float(cl["min_query_cov"])
    min_id = float(cl["min_identity"])
    min_bp = int(cl["min_cds_overlap_bp"])
    min_frac = float(cl["min_cds_overlap_frac"])

    records = harvest["records"]
    counts = quantify_flnc(cfg, harvest["transcripts_fa"]) if cl.get("quantify_flnc", True) else {}

    rows: dict[str, dict] = {}
    for tid, rec in records.items():
        rows[tid] = {
            "transcript_id": tid, "length": rec["length"],
            "cds_start": rec["cds_start"], "cds_end": rec["cds_end"],
            "cds_strand": rec["cds_strand"], "aa_len": rec["aa_len"],
            "completeness": rec["completeness"], "locus": rec["locus"],
            "alt_num": rec["alt_num"], "evg_class": rec["evg_class"], "oids": rec["oids"],
            "n_flnc_reads": counts.get(tid, 0),
            "from_reference_cds": any(o.startswith(cfg.d["tracks"]["ref_cds_prefix"])
                                      for o in str(rec["oids"]).split(",")),
            "chrom": "", "g_start": "", "g_end": "", "g_strand": "",
            "n_genome_exons": "", "query_cov": "", "identity": "",
            "cds_overlap_bp": 0, "ref_gene": "", "ref_genes_all": "",
            "antisense_bp": 0, "antisense_to_gene": "",
            "class": "novel_unmapped", "ref_evidence": "",
        }

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    seen = set()
    for rec in bam.fetch(until_eof=True):
        if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
            continue
        tid = rec.query_name
        if tid not in rows or tid in seen:
            continue
        seen.add(tid)
        row = rows[tid]
        aq, qlen, ident = aln_metrics(rec)
        qlen = qlen or row["length"]
        cov = aq / qlen if qlen else 0.0
        row.update({
            "chrom": rec.reference_name,
            "g_start": rec.reference_start + 1,
            "g_end": rec.reference_end,
            "n_genome_exons": len(aligned_segments(rec)),
            "query_cov": round(cov, 4), "identity": round(ident, 4),
        })
        aln_rev = rec.is_reverse
        row["g_strand"] = "-" if (aln_rev != (row["cds_strand"] == "-")) else "+"

        if cov < min_cov or ident < min_id:
            row["class"] = "novel_unmapped"
            row["ref_evidence"] = f"low_aln(cov={cov:.2f},id={ident:.3f})"
            continue

        # transcript CDS interval -> stored-read orientation -> genome
        cs, ce = int(row["cds_start"]) - 1, int(row["cds_end"])       # 0-based half-open
        if aln_rev:
            cs, ce = qlen - ce, qlen - cs
        blocks = project_interval(rec, cs, ce)
        cds_bp = sum(b - a for a, b in blocks)
        gstrand = row["g_strand"]

        sense = defaultdict(int)
        anti = defaultdict(int)
        for a, b in blocks:
            for s, e, strand, gid in ra.cds_overlaps(rec.reference_name, a + 1, b):
                ov = overlap_bp(a, b, s - 1, e)
                if ov <= 0:
                    continue
                (sense if strand == gstrand else anti)[gid] += ov

        best_gene, best_bp = "", 0
        if sense:
            best_gene, best_bp = max(sense.items(), key=lambda kv: kv[1])
        anti_gene, anti_bp = "", 0
        if anti:
            anti_gene, anti_bp = max(anti.items(), key=lambda kv: kv[1])
        row["antisense_to_gene"] = anti_gene
        row["cds_overlap_bp"] = best_bp
        row["antisense_bp"] = anti_bp
        row["ref_genes_all"] = ";".join(f"{g}:{b}" for g, b in
                                        sorted(sense.items(), key=lambda kv: -kv[1]))

        frac_ok = (best_bp / cds_bp >= min_frac) if cds_bp else False
        if best_bp >= min_bp and frac_ok:
            row["class"] = REF_CLASS
            row["ref_gene"] = best_gene
            row["ref_evidence"] = "cds_overlap"
        elif anti_bp >= min_bp:
            row["class"] = "novel_antisense"
        else:
            in_gene = ra.gene_overlaps(rec.reference_name,
                                       rec.reference_start + 1, rec.reference_end)
            row["class"] = "novel_intronic" if in_gene else "novel_intergenic"

    bam.close()

    # ---- reference CDS evidence carried through the evigene locus -----------
    if cl.get("propagate_ref_by_locus", True):
        locus_gene: dict[str, str] = {}
        for row in rows.values():
            for oid in str(row["oids"]).split(","):
                gid = refcds2gene.get(oid.strip())
                if gid and row["locus"]:
                    locus_gene.setdefault(row["locus"], gid)
            if row["class"] == REF_CLASS and row["locus"]:
                locus_gene.setdefault(row["locus"], row["ref_gene"])
        n_fixed = 0
        for row in rows.values():
            gid = locus_gene.get(row["locus"])
            if not gid:
                continue
            # never overrule direct genomic evidence: an antisense transcript stays
            # antisense even when evigene puts it in a reference gene's locus
            contradicted = (row["class"] == "novel_antisense" or int(row["antisense_bp"] or 0) > 0)
            if row["class"] != REF_CLASS and not contradicted:
                row["class"] = REF_CLASS
                row["ref_gene"] = gid
                row["ref_evidence"] = "evigene_locus"
                n_fixed += 1
            elif not row["ref_gene"]:
                row["ref_gene"] = gid
        if n_fixed:
            LOG.info("reclassified %d transcripts as reference via evigene locus sharing", n_fixed)

    for row in rows.values():
        for oid in str(row["oids"]).split(","):
            if oid.strip() in refcds2gene:
                row["ref_evidence"] = "reference_cds_sequence"
                row["ref_gene"] = row["ref_gene"] or refcds2gene[oid.strip()]
                row["class"] = REF_CLASS

    # ---- drop reference-CDS placeholders at loci that FLNC already supports --
    if cfg.d["tracks"].get("drop_reference_cds_when_supported", True):
        supported = {r["locus"] for r in rows.values()
                     if r["locus"] and not r["from_reference_cds"] and r["n_flnc_reads"] > 0}
        drop = [t for t, r in rows.items()
                if r["from_reference_cds"] and r["locus"] in supported]
        for t in drop:
            del rows[t]
        if drop:
            LOG.info("dropped %d reference-CDS placeholder transcripts at FLNC-supported loci",
                     len(drop))

    out = cfg.dir("05_classify")
    tsv = out / f"{cfg.sample}.annotation_table.tsv"
    cols = ["transcript_id", "class", "ref_gene", "ref_evidence", "locus", "alt_num",
            "evg_class", "length", "cds_start", "cds_end", "cds_strand", "aa_len",
            "completeness", "n_flnc_reads", "from_reference_cds", "chrom", "g_start",
            "g_end", "g_strand", "n_genome_exons", "query_cov", "identity",
            "cds_overlap_bp", "antisense_bp", "antisense_to_gene", "ref_genes_all", "oids"]
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for tid in sorted(rows):
            w.writerow(rows[tid])
    LOG.info("annotation table: %s", tsv)
    return {"table": tsv, "rows": rows, "refannot": ra}
