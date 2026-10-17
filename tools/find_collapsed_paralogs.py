#!/usr/bin/env python3
"""Find transcriptome loci that look like paralogs the reference does not resolve.

  python tools/find_collapsed_paralogs.py --config config/Em2.yaml
  python tools/find_collapsed_paralogs.py --config config/Em2.yaml --min-reads 5 --include-novel

WHY THIS IS NOT JUST "TWO LOCI ON ONE REFERENCE GENE"

evigene groups transcripts into loci from high-identity local CDS alignments, and a
paralog at ~93% CDS identity still shares plenty of >=98% stretches. On the synthetic
test set a deliberately planted paralog was placed *inside* its parent's locus as an
`alt`, not given a locus of its own. So locus membership alone cannot separate paralogs
from isoforms, and this script compares transcripts, not loci.

THE DISCRIMINATOR

Alternate isoforms of one gene differ by which exons they carry: where they do align,
they are near-identical, so protein identity stays ~100% while alignment coverage drops.
Paralogs differ by substitutions spread along the whole sequence: coverage stays high
while identity falls. The two measurements together separate the cases:

  collapsed_paralog    high coverage, identity below --same-gene-identity
                       -> two genes sharing one reference locus (assembly collapse,
                          or a family the annotation represents with one model)
  isoform_or_fragment  identity at or above the cutoff -> one gene, evigene alternates
  tandem_or_split      CDS footprints barely overlap -> distinct copies under one gene
                       model, or one model spanning two genes
  partial_overlap      alignment covers too little to judge

Verdicts are a heuristic over two numbers, not a phylogenetic test. Use the output as a
ranked list to inspect. Read support matters: with FLNC input a pair of singletons that
differ by 5% is far more likely two error-laden reads than two genes.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pysam  # noqa: E402

from egpipe.config import load_config  # noqa: E402
from egpipe.stages.s05_classify import project_interval  # noqa: E402
from egpipe.utils import LOG, iter_fasta, setup_logging  # noqa: E402


# ------------------------------------------------------------------ geometry
def merge(iv):
    if not iv:
        return []
    out = [list(x) for x in sorted(iv)]
    res = [out[0]]
    for s, e in out[1:]:
        if s <= res[-1][1]:
            res[-1][1] = max(res[-1][1], e)
        else:
            res.append([s, e])
    return [(s, e) for s, e in res]


def total(iv):
    return sum(e - s for s, e in iv)


def inter(a, b):
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if e > s:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def jaccard(a, b):
    ti = total(inter(a, b))
    tu = total(a) + total(b) - ti
    return ti / tu if tu else 0.0


# -------------------------------------------------------------------- inputs
def load_candidates(table: Path, min_reads: int):
    """Transcript rows worth comparing, keyed by transcript id."""
    keep = {}
    with open(table) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                reads = int(r["n_flnc_reads"] or 0)
            except ValueError:
                reads = 0
            if reads < min_reads:
                continue
            if not r["chrom"]:
                continue                      # unmapped: no footprint to compare
            r["_reads"] = reads
            keep[r["transcript_id"]] = r
    return keep


def cds_footprints(bam_path: Path, cand):
    out = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            r = cand.get(rec.query_name)
            if r is None or rec.query_name in out:
                continue
            qlen = rec.infer_read_length() or int(r["length"])
            cs, ce = int(r["cds_start"]) - 1, int(r["cds_end"])
            if rec.is_reverse:
                cs, ce = qlen - ce, qlen - cs
            blocks = merge(project_interval(rec, cs, ce))
            if not blocks:
                continue
            out[rec.query_name] = {
                "chrom": rec.reference_name, "blocks": blocks,
                "strand": r["g_strand"] or ("-" if rec.is_reverse else "+"),
            }
    return out


# ------------------------------------------------------------------ grouping
def build_groups(cand, foot, include_novel: bool, max_per_group: int):
    """(group_id, group_type, [transcript ids]) for anything worth comparing."""
    by_ref = defaultdict(list)
    for tid, r in cand.items():
        if tid in foot and r["class"] == "reference" and r["ref_gene"]:
            by_ref[r["ref_gene"]].append(tid)

    groups = [(g, "ref_gene", v) for g, v in by_ref.items() if len(v) > 1]

    if include_novel:
        novel = [(t, foot[t]) for t, r in cand.items()
                 if t in foot and r["class"].startswith("novel")]
        buckets = defaultdict(list)
        for t, f in novel:
            buckets[(f["chrom"], f["strand"])].append((t, f))
        n = 0
        for (chrom, strand), members in buckets.items():
            members.sort(key=lambda x: x[1]["blocks"][0][0])
            cur, cur_end = [], -1
            for t, f in members:
                s, e = f["blocks"][0][0], f["blocks"][-1][1]
                if cur and s < cur_end and total(inter(f["blocks"], cur[-1][1]["blocks"])) > 0:
                    cur.append((t, f))
                    cur_end = max(cur_end, e)
                else:
                    if len(cur) > 1:
                        n += 1
                        groups.append((f"region{n:04d}:{chrom}{strand}", "region",
                                       [c[0] for c in cur]))
                    cur, cur_end = [(t, f)], e
            if len(cur) > 1:
                n += 1
                groups.append((f"region{n:04d}:{chrom}{strand}", "region",
                               [c[0] for c in cur]))

    # cap each group by read support so a locus with hundreds of alts stays tractable
    capped = []
    for gid, gtype, tids in groups:
        if len(tids) > max_per_group:
            tids = sorted(tids, key=lambda t: -cand[t]["_reads"])[:max_per_group]
        capped.append((gid, gtype, sorted(tids)))
    return capped


# ------------------------------------------------- pairwise protein identity
def protein_identity(prot_fa: Path, wanted: set, workdir: Path, tmpdir: Path | None = None):
    if shutil.which("blastp") is None or shutil.which("makeblastdb") is None:
        LOG.warning("blastp not on PATH - cannot measure protein identity")
        return {}
    # keep blast scratch off a diskless node's RAM-backed /tmp
    import os
    env = {**os.environ}
    if tmpdir:
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
        env.update({"TMPDIR": str(tmpdir), "TMP": str(tmpdir), "TEMP": str(tmpdir)})
    sub = workdir / "paralog_candidates.faa"
    n = 0
    with open(sub, "w") as fh:
        for name, _d, seq in iter_fasta(prot_fa):
            base = name[:-6] if name.endswith("utrorf") else name
            if name in wanted or base in wanted:
                fh.write(f">{name}\n{seq}\n")
                n += 1
    if n < 2:
        return {}
    LOG.info("blastp all-vs-all over %d candidate proteins", n)
    subprocess.run(["makeblastdb", "-in", str(sub), "-dbtype", "prot",
                    "-out", str(workdir / "cand")], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    tab = workdir / "paralog_candidates.blastp.tsv"
    subprocess.run(["blastp", "-query", str(sub), "-db", str(workdir / "cand"),
                    "-outfmt", "6 qseqid sseqid pident length qlen slen bitscore",
                    "-evalue", "1e-5", "-max_target_seqs", "200",
                    "-out", str(tab)], check=True, stderr=subprocess.DEVNULL, env=env)
    best = {}
    with open(tab) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[0] == f[1]:
                continue
            q, s, pid, alen, qlen, slen = f[0], f[1], float(f[2]), int(f[3]), int(f[4]), int(f[5])
            cov = alen / min(qlen, slen) if min(qlen, slen) else 0.0
            k = tuple(sorted((q, s)))
            score = float(f[6])
            if k not in best or score > best[k][2]:
                best[k] = (pid, min(cov, 1.0), score)
    return best


# -------------------------------------------------------------------- verdict
def classify_pair(jac, pid, cov, args):
    if pid is None:
        return "no_protein_alignment"
    if jac < args.tandem_jaccard:
        return "tandem_or_split"
    if cov < args.min_aln_cov:
        return "partial_overlap"
    if pid >= args.same_gene_identity:
        return "isoform_or_fragment"
    return "collapsed_paralog"


RANK = {"collapsed_paralog": 4, "tandem_or_split": 3, "partial_overlap": 2,
        "no_protein_alignment": 1, "isoform_or_fragment": 0}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--min-reads", type=int, default=2,
                    help="skip transcripts with fewer FLNC reads (default 2)")
    ap.add_argument("--same-gene-identity", type=float, default=98.0,
                    help="protein identity >= this is one gene (default 98)")
    ap.add_argument("--min-aln-cov", type=float, default=0.70,
                    help="protein alignment must cover this much of the shorter "
                         "sequence to judge identity (default 0.70)")
    ap.add_argument("--tandem-jaccard", type=float, default=0.10,
                    help="CDS-footprint Jaccard below this is tandem/split (default 0.10)")
    ap.add_argument("--max-per-group", type=int, default=25,
                    help="compare at most N transcripts per group, by read support (25)")
    ap.add_argument("--include-novel", action="store_true",
                    help="also cluster novel-class loci by genomic region")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    setup_logging(None)
    cfg = load_config(args.config)
    table = cfg.outdir / "06_final" / f"{cfg.sample}.annotation_table.tsv"
    bam = cfg.outdir / "04_mapping" / f"{cfg.sample}.transcripts.genome.bam"
    prot = cfg.outdir / "06_final" / f"{cfg.sample}.protein.fa"
    for p in (table, bam):
        if not p.exists():
            LOG.error("missing %s - run stages 04-06 first", p)
            return 1
    out = Path(args.outdir) if args.outdir else cfg.outdir / "07_paralogs"
    out.mkdir(parents=True, exist_ok=True)

    cand = load_candidates(table, args.min_reads)
    LOG.info("%d transcripts pass --min-reads %d and are mapped", len(cand), args.min_reads)
    foot = cds_footprints(bam, cand)
    LOG.info("CDS footprints for %d transcripts", len(foot))
    groups = build_groups(cand, foot, args.include_novel, args.max_per_group)
    LOG.info("%d groups with more than one transcript", len(groups))
    if not groups:
        print("\nno reference gene or region carries more than one transcript")
        return 0

    wanted = {t for _g, _ty, ts in groups for t in ts}
    # blast intermediates go in a subdirectory so the results directory holds results
    blastdir = out / "blast"
    blastdir.mkdir(parents=True, exist_ok=True)
    pid_map = (protein_identity(prot, wanted, blastdir, cfg.tmpdir)
               if prot.exists() else {})

    pair_rows, group_rows = [], []
    for gid, gtype, tids in groups:
        verdicts = []
        for i in range(len(tids)):
            for j in range(i + 1, len(tids)):
                a, b = tids[i], tids[j]
                fa, fb = foot[a], foot[b]
                jac = jaccard(fa["blocks"], fb["blocks"]) if fa["chrom"] == fb["chrom"] else 0.0
                hit = pid_map.get(tuple(sorted((a, b))))
                pid = hit[0] if hit else None
                cov = hit[1] if hit else 0.0
                v = classify_pair(jac, pid, cov, args)
                verdicts.append(v)
                pair_rows.append({
                    "group": gid, "group_type": gtype,
                    "transcript_a": a, "transcript_b": b,
                    "locus_a": cand[a]["locus"], "locus_b": cand[b]["locus"],
                    "same_evigene_locus": "yes" if cand[a]["locus"] == cand[b]["locus"] else "no",
                    "reads_a": cand[a]["_reads"], "reads_b": cand[b]["_reads"],
                    "aa_a": cand[a]["aa_len"], "aa_b": cand[b]["aa_len"],
                    "protein_identity": "" if pid is None else round(pid, 2),
                    "protein_aln_cov": round(cov, 3),
                    "cds_footprint_jaccard": round(jac, 4),
                    "genome_identity_a": cand[a]["identity"],
                    "genome_identity_b": cand[b]["identity"],
                    "span_a": f"{fa['chrom']}:{fa['blocks'][0][0]+1}-{fa['blocks'][-1][1]}{fa['strand']}",
                    "span_b": f"{fb['chrom']}:{fb['blocks'][0][0]+1}-{fb['blocks'][-1][1]}{fb['strand']}",
                    "verdict": v,
                })
        top = max(verdicts, key=lambda v: RANK[v]) if verdicts else ""
        group_rows.append({
            "group": gid, "group_type": gtype, "n_transcripts": len(tids),
            "n_evigene_loci": len({cand[t]["locus"] for t in tids}),
            "total_reads": sum(cand[t]["_reads"] for t in tids),
            "group_verdict": top, "n_pairs": len(verdicts),
            "n_collapsed_paralog": verdicts.count("collapsed_paralog"),
            "n_tandem_or_split": verdicts.count("tandem_or_split"),
            "n_isoform_or_fragment": verdicts.count("isoform_or_fragment"),
            "min_protein_identity": min(
                [p["protein_identity"] for p in pair_rows
                 if p["group"] == gid and p["protein_identity"] != ""], default=""),
        })

    group_rows.sort(key=lambda r: (-RANK.get(r["group_verdict"], 0), -r["total_reads"]))
    gp = out / f"{cfg.sample}.paralog_groups.tsv"
    pp = out / f"{cfg.sample}.paralog_pairs.tsv"
    hist_txt = out / f"{cfg.sample}.identity_histogram.txt"
    hist_tsv = out / f"{cfg.sample}.identity_histogram.tsv"
    # Missing numerics are written as NA, not "", so read.delim() in R types the
    # column numeric instead of character.
    for path, data in ((gp, group_rows), (pp, pair_rows)):
        with open(path, "w", newline="") as fh:
            if data:
                w = csv.DictWriter(fh, fieldnames=list(data[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows([{k: ("NA" if v == "" else v) for k, v in r.items()}
                             for r in data])

    # The identity threshold is the whole call, so show where the pairs actually sit
    # rather than asking anyone to trust the default. Alternates driven by read error
    # pile up just under 100; a separate mode lower down is the paralog signal. Bars are
    # log-scaled because the isoform mode is usually orders of magnitude taller than the
    # paralog mode, which would otherwise be invisible.
    pids = sorted(p["protein_identity"] for p in pair_rows if p["protein_identity"] != "")
    if pids:
        import math
        lo = max(50.0, math.floor(min(pids)))
        step = 1.0 if (100 - lo) <= 45 else 2.0
        bins = []
        x = lo
        while x < 100:
            hi = min(100.0, x + step)
            bins.append((x, hi, sum(1 for v in pids if x <= v < hi or (hi >= 100 and v >= 100))))
            x += step
        peak = max((c for _a, _b, c in bins), default=0) or 1
        q = lambda f: pids[min(len(pids) - 1, int(len(pids) * f))]
        lines = [f"protein identity of compared pairs (n={len(pids)}, "
                 f"--same-gene-identity={args.same_gene_identity}, bars log-scaled)"]
        for a, b, c in bins:
            # empty bins are kept: the trough between the isoform mode and any paralog
            # mode is the thing you are looking for, and suppressing zeros hides it
            width = int(40 * math.log10(c + 1) / math.log10(peak + 1)) if c else 0
            mark = "  <-- current cutoff" if a <= args.same_gene_identity < b else ""
            lines.append(f"  {a:5.1f}-{b:<5.1f} {'#' * width:<40} {c}{mark}")
        lines.append(f"  min {pids[0]:.2f}   p1 {q(0.01):.2f}   p5 {q(0.05):.2f}   "
                     f"p50 {q(0.50):.2f}   max {pids[-1]:.2f}")
        print()
        print("\n".join(lines))

        hist_txt.write_text(
            "\n".join(lines)
            + "\n\n"
            + f"# bin counts: {hist_tsv.name}\n"
            + f"# per-pair values for re-binning: {pp.name} (protein_identity column)\n"
            + f"# min_reads={args.min_reads} min_aln_cov={args.min_aln_cov} "
              f"tandem_jaccard={args.tandem_jaccard} max_per_group={args.max_per_group}\n"
        )
        with open(hist_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["bin_low", "bin_high", "count"])
            for a, b, c in bins:
                w.writerow([f"{a:.1f}", f"{b:.1f}", c])
        print(f"  histogram  -> {hist_txt.name}")
        print(f"  bin counts -> {hist_tsv.name}")
        print(f"  raw values -> {pp.name} (protein_identity, protein_aln_cov, "
              f"cds_footprint_jaccard)")

    counts = defaultdict(int)
    for r in group_rows:
        counts[r["group_verdict"]] += 1
    print(f"\ngroups examined: {len(group_rows)}   pairs: {len(pair_rows)}")
    for k in sorted(RANK, key=lambda x: -RANK[x]):
        print(f"  {k:<22} {counts.get(k, 0)}")
    print(f"\nwrote {out}/")
    top = [r for r in group_rows if r["group_verdict"] == "collapsed_paralog"]
    if top:
        print("\ncollapsed-paralog candidates (highest read support first):")
        for r in top[:15]:
            print(f"  {r['group']:<24} transcripts={r['n_transcripts']} "
                  f"loci={r['n_evigene_loci']} reads={r['total_reads']:<6} "
                  f"min_prot_id={r['min_protein_identity']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
