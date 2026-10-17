"""Stage 06 - write the final GFF3 in TRANSCRIPT coordinates, plus deliverables.

Each transcript is its own GFF sequence region: seqid = transcript ID, one exon
spanning the whole transcript, and the CDS placed at the evigene ORF offsets.
A genome-coordinate GFF is written alongside as a cross-check / browsing aid.
"""
from __future__ import annotations

import csv
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import pysam

from ..utils import LOG, iter_fasta, write_fasta
from .s05_classify import aligned_segments, project_interval

STAGE = "06_annotate"


def _attr(d: dict) -> str:
    parts = []
    for k, v in d.items():
        if v is None or v == "":
            continue
        parts.append(f"{k}={str(v).replace(';', '%3B').replace('=', '%3D').replace(',', '%2C')}")
    return ";".join(parts)


def write_transcript_gff(cfg, rows: dict, path: Path) -> None:
    src = cfg.d["output"]["gff_source"]
    emit_utr = bool(cfg.d["output"]["emit_utr"])
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"# transcriptome coordinates ({cfg.sample}); seqid = transcript ID\n")
        for tid in sorted(rows):
            r = rows[tid]
            ln = int(r["length"])
            cs, ce = int(r["cds_start"]), int(r["cds_end"])
            strand = r["cds_strand"] or "+"
            fh.write(f"##sequence-region {tid} 1 {ln}\n")
            common = {
                "class": r["class"],
                "ref_gene": r["ref_gene"],
                "ref_evidence": r["ref_evidence"],
                "evigene_locus": r["locus"],
                "evigene_class": r["evg_class"],
                "alt_num": r["alt_num"],
                "aa_len": r["aa_len"],
                "completeness": r["completeness"],
                "n_flnc_reads": r["n_flnc_reads"],
                "from_reference_cds": "true" if r["from_reference_cds"] else None,
            }
            if r["chrom"]:
                common["genome_locus"] = f"{r['chrom']}:{r['g_start']}-{r['g_end']}{r['g_strand']}"
                common["genome_identity"] = r["identity"]
                common["genome_query_cov"] = r["query_cov"]
                common["genome_exons"] = r["n_genome_exons"]
            gene_id = f"{tid}.gene"
            fh.write("\t".join([tid, src, "gene", "1", str(ln), ".", strand, ".",
                                _attr({"ID": gene_id, **common})]) + "\n")
            fh.write("\t".join([tid, src, "mRNA", "1", str(ln), ".", strand, ".",
                                _attr({"ID": tid, "Parent": gene_id, **common})]) + "\n")
            fh.write("\t".join([tid, src, "exon", "1", str(ln), ".", strand, ".",
                                _attr({"ID": f"{tid}.exon1", "Parent": tid})]) + "\n")
            if emit_utr and cs > 1:
                feat = "five_prime_UTR" if strand == "+" else "three_prime_UTR"
                fh.write("\t".join([tid, src, feat, "1", str(cs - 1), ".", strand, ".",
                                    _attr({"ID": f"{tid}.{feat}", "Parent": tid})]) + "\n")
            fh.write("\t".join([tid, src, "CDS", str(cs), str(ce), ".", strand, "0",
                                _attr({"ID": f"{tid}.cds", "Parent": tid,
                                       "completeness": r["completeness"]})]) + "\n")
            if emit_utr and ce < ln:
                feat = "three_prime_UTR" if strand == "+" else "five_prime_UTR"
                fh.write("\t".join([tid, src, feat, str(ce + 1), str(ln), ".", strand, ".",
                                    _attr({"ID": f"{tid}.{feat}", "Parent": tid})]) + "\n")
    LOG.info("transcript-coordinate GFF: %s", path)


def write_genome_gff(cfg, rows: dict, bam_path: Path, path: Path) -> None:
    src = cfg.d["output"]["gff_source"]
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    seen = set()
    n = 0
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"# genome coordinates ({cfg.sample}); cross-check for the transcript GFF\n")
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            tid = rec.query_name
            if tid not in rows or tid in seen:
                continue
            seen.add(tid)
            r = rows[tid]
            strand = r["g_strand"] or ("-" if rec.is_reverse else "+")
            blocks = aligned_segments(rec)
            if not blocks:
                continue
            qlen = rec.infer_read_length() or int(r["length"])
            cs, ce = int(r["cds_start"]) - 1, int(r["cds_end"])
            if rec.is_reverse:
                cs, ce = qlen - ce, qlen - cs
            cds_blocks = project_interval(rec, cs, ce)
            gid = f"{tid}.gene"
            attrs = {"ID": gid, "class": r["class"], "ref_gene": r["ref_gene"],
                     "evigene_locus": r["locus"]}
            fh.write("\t".join([rec.reference_name, src, "gene",
                                str(rec.reference_start + 1), str(rec.reference_end),
                                ".", strand, ".", _attr(attrs)]) + "\n")
            fh.write("\t".join([rec.reference_name, src, "mRNA",
                                str(rec.reference_start + 1), str(rec.reference_end),
                                ".", strand, ".",
                                _attr({"ID": tid, "Parent": gid, **{k: v for k, v in attrs.items() if k != "ID"}})]) + "\n")
            for i, (a, b, ra_, rb) in enumerate(blocks, 1):
                fh.write("\t".join([rec.reference_name, src, "exon", str(ra_ + 1), str(rb),
                                    ".", strand, ".",
                                    _attr({"ID": f"{tid}.e{i}", "Parent": tid})]) + "\n")
            for i, (a, b) in enumerate(cds_blocks, 1):
                fh.write("\t".join([rec.reference_name, src, "CDS", str(a + 1), str(b),
                                    ".", strand, ".",
                                    _attr({"ID": f"{tid}.c{i}", "Parent": tid})]) + "\n")
            n += 1
    bam.close()
    LOG.info("genome-coordinate GFF: %s (%d mapped transcripts)", path, n)


def write_summary(cfg, rows: dict, ra, path: Path) -> None:
    cls = Counter(r["class"] for r in rows.values())
    loci_by_class = defaultdict(set)
    for r in rows.values():
        loci_by_class[r["class"]].add(r["locus"] or r["transcript_id"])
    ref_genes_hit = {r["ref_gene"] for r in rows.values() if r["ref_gene"]}
    all_ref_genes = {g for chrom in ra.genes for (_s, _e, _st, g) in ra.genes[chrom]}
    novel_loci = {r["locus"] or r["transcript_id"] for r in rows.values()
                  if r["class"].startswith("novel")}
    # reference genes represented ONLY by the reference-CDS track (no FLNC support)
    supported_genes = {r["ref_gene"] for r in rows.values()
                       if r["ref_gene"] and not r["from_reference_cds"]}
    refcds_only_genes = {r["ref_gene"] for r in rows.values()
                         if r["ref_gene"] and r["from_reference_cds"]} - supported_genes

    lines = [
        f"# {cfg.sample} EvidentialGene transcriptome summary",
        "",
        f"transcripts_total\t{len(rows)}",
        f"evigene_loci_total\t{len({r['locus'] or r['transcript_id'] for r in rows.values()})}",
        "",
        "## transcripts by class",
    ]
    for k in ("reference", "novel_intergenic", "novel_intronic", "novel_antisense",
              "novel_unmapped"):
        lines.append(f"{k}\t{cls.get(k, 0)}\t(loci: {len(loci_by_class.get(k, ()))})")
    lines += [
        "",
        "## reference recovery",
        f"reference_genes_in_annotation\t{len(all_ref_genes)}",
        f"reference_genes_recovered\t{len(ref_genes_hit & all_ref_genes)}",
        f"reference_genes_missing\t{len(all_ref_genes - ref_genes_hit)}",
        f"reference_genes_without_flnc_support\t{len(refcds_only_genes)}"
        "\t(appended from the annotation, no FLNC-supported transcript)",
        "",
        "## novelty",
        f"novel_transcripts\t{sum(v for k, v in cls.items() if k.startswith('novel'))}",
        f"novel_loci\t{len(novel_loci)}",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")
    LOG.info("summary:\n%s", "\n".join(lines))

    missing = sorted(all_ref_genes - ref_genes_hit)
    (path.parent / f"{cfg.sample}.reference_genes_missing.txt").write_text(
        "\n".join(missing) + ("\n" if missing else ""))


def append_unrecovered_reference(cfg, rows: dict, ra) -> dict:
    """Add reference genes the reads never recovered, as CDS-only transcripts.

    Keeping the reference out of tr2aacds and adding it here avoids a real failure
    mode: after cluster2 the FLNC consensus CDS is often byte-identical to the
    annotated CDS, evigene's perfect-duplicate step then prefers the CDS-only
    sequence (100% coding), and the full-length transcript with its UTRs is dropped.
    """
    if not cfg.d["transcriptome"].get("add_unrecovered_reference", True):
        return {}
    inp = cfg.outdir / "01_inputs"
    cds_fa, prot_fa = inp / "ref_cds.raw.fa", inp / "ref_prot.raw.fa"
    if not cds_fa.exists():
        LOG.warning("no %s - cannot append unrecovered reference genes", cds_fa)
        return {}

    recovered = {r["ref_gene"] for r in rows.values()
                 if r["class"] == "reference" and r["ref_gene"]}
    span = {}
    for chrom, items in ra.genes.items():
        for gs, ge, strand, gid in items:
            span[gid] = (chrom, gs, ge, strand)

    best: dict[str, tuple[str, str]] = {}
    for name, _desc, seq in iter_fasta(cds_fa):
        gid = ra.tx2gene.get(name, name)
        if gid in recovered or not seq:
            continue
        if gid not in best or len(seq) > len(best[gid][1]):
            best[gid] = (name, seq)
    if not best:
        LOG.info("every reference gene was recovered from the reads")
        return {}

    prot = {}
    if prot_fa.exists():
        wanted = {tid for tid, _s in best.values()}
        prot = {n: sq for n, _d, sq in iter_fasta(prot_fa) if n in wanted}

    extra = {}
    for gid, (tid, seq) in sorted(best.items()):
        chrom, gs, ge, strand = span.get(gid, ("", "", "", "+"))
        rows[tid] = {
            "transcript_id": tid, "length": len(seq),
            "cds_start": 1, "cds_end": len(seq), "cds_strand": "+",
            "aa_len": max(0, len(seq) // 3 - 1), "completeness": "reference_cds",
            "locus": gid, "alt_num": 1, "evg_class": "reference_only",
            "oids": "", "n_flnc_reads": 0, "from_reference_cds": True,
            "chrom": chrom, "g_start": gs, "g_end": ge, "g_strand": strand,
            "n_genome_exons": "", "query_cov": "", "identity": "",
            "cds_overlap_bp": len(seq), "ref_gene": gid,
            "ref_genes_all": f"{gid}:{len(seq)}", "antisense_bp": 0,
            "antisense_to_gene": "", "class": "reference",
            "ref_evidence": "reference_only_no_flnc",
        }
        extra[tid] = (seq, prot.get(tid, ""))
    LOG.info("appended %d reference genes with no FLNC support", len(extra))
    return extra


def run_stage(cfg, harvest: dict, classified: dict, bam_path: Path) -> Path:
    rows = classified["rows"]
    ra = classified["refannot"]
    extra = append_unrecovered_reference(cfg, rows, ra)
    final = cfg.dir("06_final")

    gff = final / f"{cfg.sample}.transcriptome.gff3"
    write_transcript_gff(cfg, rows, gff)

    if cfg.d["classify"].get("emit_genome_gff", True):
        write_genome_gff(cfg, rows, bam_path, final / f"{cfg.sample}.genome.gff3")

    # subset the sequence files to the transcripts that survived classification
    for src_name, dest_name in (
        (f"{cfg.sample}.transcripts.fa", f"{cfg.sample}.transcriptome.fa"),
        (f"{cfg.sample}.cds.fa", f"{cfg.sample}.cds.fa"),
        (f"{cfg.sample}.protein.fa", f"{cfg.sample}.protein.fa"),
    ):
        src = cfg.outdir / "03_transcriptome" / src_name
        if not src.exists():
            continue
        n = 0
        with open(final / dest_name, "w") as fh:
            for name, desc, seq in iter_fasta(src):
                base = name[:-6] if name.endswith("utrorf") else name
                if name in rows or base in rows:
                    write_fasta(fh, name, seq, desc)
                    n += 1
            for tid, (cds_seq, prot_seq) in extra.items():
                if dest_name.endswith("protein.fa"):
                    if prot_seq:
                        write_fasta(fh, tid, prot_seq,
                                    f"type=protein; source=reference_annotation; gene={rows[tid]['locus']}")
                        n += 1
                else:
                    write_fasta(fh, tid, cds_seq,
                                f"type=cds; source=reference_annotation; "
                                f"gene={rows[tid]['locus']}; clen={len(cds_seq)}; "
                                f"strand=+; offs=1-{len(cds_seq)};")
                    n += 1
        LOG.info("wrote %s (%d seqs)", dest_name, n)

    # rewrite the table so it includes the appended reference genes
    import csv as _csv
    cols = ["transcript_id", "class", "ref_gene", "ref_evidence", "locus", "alt_num",
            "evg_class", "length", "cds_start", "cds_end", "cds_strand", "aa_len",
            "completeness", "n_flnc_reads", "from_reference_cds", "chrom", "g_start",
            "g_end", "g_strand", "n_genome_exons", "query_cov", "identity",
            "cds_overlap_bp", "antisense_bp", "antisense_to_gene", "ref_genes_all", "oids"]
    with open(final / f"{cfg.sample}.annotation_table.tsv", "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for tid in sorted(rows):
            w.writerow(rows[tid])
    write_summary(cfg, rows, ra, final / f"{cfg.sample}.summary.txt")
    LOG.info("final outputs in %s", final)
    return final
