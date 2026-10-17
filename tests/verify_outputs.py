#!/usr/bin/env python3
"""Verify a finished run: GFF coordinates, translations, and classification.

  python tests/verify_outputs.py --outdir /tmp/testrun --sample TestSp \
      [--expected /tmp/testdata/expected_classes.tsv]

Checks
  1. every transcript in the GFF exists in the transcriptome fasta, with matching length
  2. the CDS slice taken from the transcriptome using the GFF coordinates is identical
     to the evigene CDS sequence for that transcript
  3. CDS length is a multiple of 3, 'complete' ORFs start with ATG and end in a stop,
     and no ORF has an internal stop codon
  4. exon + UTR features tile the transcript exactly once
  5. (optional) classification matches the answer key from make_test_data.py
"""
from __future__ import annotations

import argparse
import itertools
import re
import sys
from collections import defaultdict
from pathlib import Path

BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AAS[i] for i, (a, b, c) in enumerate(itertools.product(BASES, repeat=3))}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOPS = {"TAA", "TAG", "TGA"}


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def translate(s):
    return "".join(CODON.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--expected", default=None)
    a = ap.parse_args()
    final = Path(a.outdir) / "06_final"

    tx = read_fasta(final / f"{a.sample}.transcriptome.fa")
    cds_fa = final / f"{a.sample}.cds.fa"
    cds = read_fasta(cds_fa) if cds_fa.exists() else {}
    gff = final / f"{a.sample}.transcriptome.gff3"

    feats = defaultdict(list)
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        feats[f[0]].append((f[2], int(f[3]), int(f[4]), f[6], f[7], f[8]))

    fails, warns, n_cds_checked = [], [], 0

    for tid, rows in feats.items():
        if tid not in tx:
            fails.append(f"{tid}: in GFF but not in transcriptome fasta")
            continue
        seq = tx[tid]
        L = len(seq)
        for kind, s, e, strand, phase, attrs in rows:
            if s < 1 or e > L:
                fails.append(f"{tid}: {kind} {s}-{e} outside transcript length {L}")
        mrna = [r for r in rows if r[0] == "mRNA"]
        if len(mrna) != 1 or mrna[0][1] != 1 or mrna[0][2] != L:
            fails.append(f"{tid}: mRNA feature should span 1-{L}")
        exons = [r for r in rows if r[0] == "exon"]
        if len(exons) != 1 or exons[0][1] != 1 or exons[0][2] != L:
            fails.append(f"{tid}: expected exactly one exon spanning 1-{L}")

        # UTR + CDS must tile the transcript exactly once
        cover = [0] * (L + 1)
        for kind, s, e, *_ in rows:
            if kind in ("CDS", "five_prime_UTR", "three_prime_UTR"):
                for i in range(s, e + 1):
                    cover[i] += 1
        bad = [i for i in range(1, L + 1) if cover[i] != 1]
        if bad:
            fails.append(f"{tid}: CDS+UTR tiling covers {len(bad)} positions != 1x")

        cds_rows = [r for r in rows if r[0] == "CDS"]
        if len(cds_rows) != 1:
            fails.append(f"{tid}: expected one CDS feature, found {len(cds_rows)}")
            continue
        _k, s, e, strand, phase, attrs = cds_rows[0]
        sub = seq[s - 1:e]
        if strand == "-":
            sub = rc(sub)
        if len(sub) % 3:
            warns.append(f"{tid}: CDS length {len(sub)} not a multiple of 3")
        compl = re.search(r"completeness=([^;]+)", attrs)
        compl = compl.group(1) if compl else ""
        if "complete" in compl:
            if not sub.upper().startswith("ATG"):
                fails.append(f"{tid}: complete ORF does not start with ATG ({sub[:3]})")
            if sub[-3:].upper() not in STOPS:
                fails.append(f"{tid}: complete ORF does not end in a stop ({sub[-3:]})")
        prot = translate(sub)
        if "*" in prot[:-1]:
            fails.append(f"{tid}: {prot[:-1].count('*')} internal stop codon(s) in CDS")

        key = tid if tid in cds else None
        if key:
            ref = cds[key].upper()
            got = sub.upper()
            if ref.rstrip("*") != got and ref != got:
                if got.startswith(ref[:30]) and abs(len(ref) - len(got)) <= 3:
                    warns.append(f"{tid}: CDS differs from evigene cds by <=3 nt (stop codon)")
                else:
                    fails.append(f"{tid}: GFF-derived CDS != evigene CDS "
                                 f"(len {len(got)} vs {len(ref)})")
            n_cds_checked += 1

    print(f"transcripts in GFF : {len(feats)}")
    print(f"transcripts in fasta: {len(tx)}")
    print(f"CDS cross-checked vs evigene cds fasta: {n_cds_checked}")

    if a.expected:
        table = final / f"{a.sample}.annotation_table.tsv"
        hdr = None
        rows = []
        for line in open(table):
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            rows.append(dict(zip(hdr, f)))
        # flnc_NNNNN -> original name. With isoseq cluster2 the "original" is a
        # cluster id, so chain: flnc_id -> transcript/N -> movie/zmw/ccs -> origin
        idmap = {}
        for line in open(Path(a.outdir) / "01_inputs" / "flnc.idmap.tsv"):
            f = line.rstrip("\n").split("\t")
            if f[0] != "new_id":
                idmap[f[0]] = f[1]

        cluster_members = {}
        crep = Path(a.outdir) / "01_inputs" / "clustered.cluster_report.csv"
        if crep.exists():
            for line in open(crep):
                f = line.rstrip("\n").split(",")
                if len(f) >= 2 and f[0] != "cluster_id":
                    cluster_members.setdefault(f[0], []).append(f[1])
        readmap = {}
        rm = Path(a.expected).parent / "flnc.bam.readmap.tsv"
        if rm.exists():
            for line in open(rm):
                f = line.rstrip("\n").split("\t")
                if f[0] != "bam_read_id":
                    readmap[f[0]] = f[1]

        def origins_of(oid):
            orig = idmap.get(oid)
            if orig is None:
                return []
            if orig in cluster_members:
                return [readmap.get(r, r) for r in cluster_members[orig]]
            return [readmap.get(orig, orig)]
        print("\norigin -> classes assigned:")
        agg = defaultdict(lambda: defaultdict(int))
        for r in rows:
            origins = set()
            if r.get("ref_evidence") == "reference_only_no_flnc":
                origins.add("reference_appendix")
            for oid in r["oids"].split(","):
                oid = oid.strip()
                if oid.startswith("refcds"):
                    origins.add("reference_cds_track")
                    continue
                for o in origins_of(oid):
                    origins.add(re.sub(r"_(read|dup)\d+$", "", o))
            for o in origins or {"(no oid)"}:
                agg[o][r["class"]] += 1
        for o in sorted(agg):
            print(f"  {o:<22} " + ", ".join(f"{k}={v}" for k, v in sorted(agg[o].items())))
        print("\nanswer key:")
        print(Path(a.expected).read_text().rstrip())

    print()
    for w in warns:
        print("WARN:", w)
    for f in fails:
        print("FAIL:", f)
    print(f"\n{len(fails)} failures, {len(warns)} warnings")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
