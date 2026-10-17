#!/usr/bin/env python3
"""Build a small synthetic dataset that exercises every classification branch.

Creates: genome.fa, reference.gff3, flnc.fa, flnc.bam (if samtools is present),
and expected_classes.tsv holding the answer key.
"""
from __future__ import annotations

import argparse
import random
import subprocess
import shutil
from pathlib import Path

CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"]
STOPS = {"TAA", "TAG", "TGA"}
SENSE = [c for c in CODONS if c not in STOPS]
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def rand_seq(rng, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_cds(rng, n_codons: int) -> str:
    return "ATG" + "".join(rng.choice(SENSE) for _ in range(n_codons - 2)) + "TAA"


def diverge_cds(rng, cds: str, codon_rate: float) -> str:
    """Diverge a CDS by swapping whole codons, so no stop codon is ever introduced.

    Swapping a codon changes ~2 of its 3 bases, so nucleotide divergence lands near
    (2/3 * codon_rate) and the protein diverges by about codon_rate.
    """
    out = [cds[i:i + 3] for i in range(0, len(cds), 3)]
    for i in range(1, len(out) - 1):          # keep ATG and the stop
        if rng.random() < codon_rate:
            out[i] = rng.choice(SENSE)
    return "".join(out)


def mutate(rng, seq: str, rate: float) -> str:
    if rate <= 0:
        return seq
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


class Gene:
    """A spliced gene laid out in transcript orientation, then placed on a strand."""

    def __init__(self, rng, name, n_codons=220, utr5=110, utr3=160, n_exons=3,
                 intron_range=(280, 460)):
        self.name = name
        cds = make_cds(rng, n_codons)
        self.mrna = rand_seq(rng, utr5) + cds + rand_seq(rng, utr3)
        self.cds_start = utr5 + 1                      # 1-based, transcript coords
        self.cds_end = utr5 + len(cds)
        cuts = sorted(rng.sample(range(180, len(self.mrna) - 180), n_exons - 1))
        bounds = [0] + cuts + [len(self.mrna)]
        self.exons_tx = [(bounds[i] + 1, bounds[i + 1]) for i in range(n_exons)]
        self.introns = [("GT" + rand_seq(rng, rng.randint(*intron_range)) + "AG")
                        for _ in range(n_exons - 1)]
        parts, self.exon_u = [], []
        pos = 0
        for i, (s, e) in enumerate(self.exons_tx):
            seg = self.mrna[s - 1:e]
            parts.append(seg)
            self.exon_u.append((pos + 1, pos + len(seg)))
            pos += len(seg)
            if i < len(self.introns):
                parts.append(self.introns[i])
                pos += len(self.introns[i])
        self.unspliced = "".join(parts)

    def place(self, chrom, offset, strand):
        """Return (genomic_seq, exon_intervals, strand) with 1-based genome coords."""
        self.chrom, self.strand = chrom, strand
        L = len(self.unspliced)
        if strand == "+":
            g = self.unspliced
            self.exons_g = [(offset + s, offset + e) for s, e in self.exon_u]
        else:
            g = rc(self.unspliced)
            self.exons_g = [(offset + (L - e) + 1, offset + (L - s) + 1)
                            for s, e in self.exon_u]
        self.exons_g.sort()
        self.offset = offset
        return g

    def tx_to_genome(self, ts, te):
        """Map a transcript interval (1-based inclusive) to genomic blocks."""
        blocks = []
        for (es, ee), (gs, ge) in zip(self.exons_tx, self.exon_u):
            lo, hi = max(es, ts), min(ee, te)
            if hi < lo:
                continue
            us = gs + (lo - es)
            ue = gs + (hi - es)
            L = len(self.unspliced)
            if self.strand == "+":
                blocks.append((self.offset + us, self.offset + ue))
            else:
                blocks.append((self.offset + (L - ue) + 1, self.offset + (L - us) + 1))
        blocks.sort()
        return blocks


def gff_lines(g: Gene):
    out = []
    gs = min(s for s, _ in g.exons_g)
    ge = max(e for _, e in g.exons_g)
    gid, tid = g.name, g.name + ".t1"
    out.append(f"{g.chrom}\ttest\tgene\t{gs}\t{ge}\t.\t{g.strand}\t.\tID={gid}")
    out.append(f"{g.chrom}\ttest\tmRNA\t{gs}\t{ge}\t.\t{g.strand}\t.\tID={tid};Parent={gid}")
    for i, (s, e) in enumerate(g.exons_g, 1):
        out.append(f"{g.chrom}\ttest\texon\t{s}\t{e}\t.\t{g.strand}\t.\tID={tid}.e{i};Parent={tid}")
    for i, (s, e) in enumerate(g.tx_to_genome(g.cds_start, g.cds_end), 1):
        out.append(f"{g.chrom}\ttest\tCDS\t{s}\t{e}\t.\t{g.strand}\t0\tID={tid}.c{i};Parent={tid}")
    return out


def write_pacbio_bam(reads, path: Path) -> None:
    """Write a PacBio-style unaligned FLNC BAM.

    The @RG with PL:PACBIO and DS:READTYPE=CCS matters: `isoseq cluster2` hangs
    with no output on a uBAM that lacks it (e.g. one made by `samtools import`).
    """
    try:
        import pysam
    except ImportError:
        print("pysam not available - skipping BAM output")
        return
    movie = "m84001_260101_000000_s1"
    rgid = "a1b2c3d4"
    header = {
        "HD": {"VN": "1.6", "SO": "unknown", "pb": "5.0.0"},
        "RG": [{"ID": rgid, "PL": "PACBIO", "PM": "REVIO", "PU": movie,
                "DS": ("READTYPE=CCS;BINDINGKIT=102-739-100;SEQUENCINGKIT=102-118-800;"
                       "BASECALLERVERSION=1.0;FRAMERATEHZ=100.000000")}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for i, (_name, seq) in enumerate(reads):
            a = pysam.AlignedSegment()
            a.query_name = f"{movie}/{1000 + i}/ccs"
            a.flag = 4
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("~" * len(seq))
            a.set_tag("np", 20, "i")
            a.set_tag("rq", 0.999, "f")
            a.set_tag("zm", 1000 + i, "i")
            a.set_tag("RG", rgid, "Z")
            fh.write(a)
    # map the synthetic read names onto the movie/zmw ids for the answer key
    with open(str(path) + ".readmap.tsv", "w") as fh:
        fh.write("bam_read_id\torigin_read\n")
        for i, (name, _seq) in enumerate(reads):
            fh.write(f"{movie}/{1000 + i}/ccs\t{name}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="testdata")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    chroms = {"chr1": [], "chr2": []}
    layout = {"chr1": 0, "chr2": 0}
    gff, reads, key = [], [], []

    def add_spacer(chrom, n):
        chroms[chrom].append(rand_seq(rng, n))
        layout[chrom] += n

    def place(gene, chrom, strand):
        add_spacer(chrom, rng.randint(1200, 2200))
        seq = gene.place(chrom, layout[chrom], strand)
        chroms[chrom].append(seq)
        layout[chrom] += len(seq)
        return gene

    def emit_reads(gene, n, label, rate=0.002, exact_dups=0):
        for i in range(n):
            reads.append((f"{gene.name}_read{i}", mutate(rng, gene.mrna, rate)))
        for i in range(exact_dups):
            reads.append((f"{gene.name}_dup{i}", gene.mrna))
        key.append((gene.name, label))

    # --- reference genes, 4 with FLNC support, 1 without --------------------
    annotated = []
    for i in range(1, 6):
        strand = "-" if i == 3 else "+"
        irange = (900, 1100) if i == 2 else (280, 460)   # gene2 hosts an intronic novel gene
        g = place(Gene(rng, f"refgene{i}", n_codons=200 + 20 * i,
                       intron_range=irange), "chr1", strand)
        annotated.append(g)
        gff += gff_lines(g)
        if i <= 4:
            emit_reads(g, 4, "reference", exact_dups=2 if i == 1 else 0)
        else:
            key.append((g.name, "reference(no FLNC; from reference CDS track)"))

    # --- novel intergenic locus on chr2 (not in the GFF) -------------------
    nov = place(Gene(rng, "novel_intergenic1", n_codons=260), "chr2", "+")
    emit_reads(nov, 5, "novel_intergenic")

    # --- novel gene inside an intron of refgene2, same strand --------------
    intronic_host = annotated[1]
    intronic = Gene(rng, "novel_intronic1", n_codons=70, utr5=30, utr3=30, n_exons=1)
    # rebuild the host so the new gene sits inside its first intron
    ih_start = intronic_host.offset + intronic_host.exon_u[0][1] + 100
    host_intron = intronic_host.introns[0]
    assert len(intronic.unspliced) < len(host_intron) - 200, "host intron too short"
    if True:
        chunk = list("".join(chroms["chr1"]))
        chunk[ih_start:ih_start + len(intronic.unspliced)] = list(intronic.unspliced)
        chroms["chr1"] = ["".join(chunk)]
        intronic.chrom, intronic.strand, intronic.offset = "chr1", "+", ih_start
        intronic.exons_g = [(ih_start + 1, ih_start + len(intronic.unspliced))]
        emit_reads(intronic, 4, "novel_intronic")

    # --- antisense transcript over refgene1's CDS --------------------------
    host = annotated[0]
    anti_mrna = rc(host.mrna)
    reads += [(f"antisense1_read{i}", mutate(rng, anti_mrna, 0.002)) for i in range(4)]
    key.append(("antisense1", "novel_antisense (reverse complement of refgene1 mRNA)"))

    # --- a transcript that is not in the genome at all ---------------------
    off = Gene(rng, "offgenome1", n_codons=180)
    reads += [(f"offgenome1_read{i}", mutate(rng, off.mrna, 0.002)) for i in range(4)]
    key.append(("offgenome1", "novel_unmapped (not present in the genome)"))

    # --- collapsed paralog: two divergent loci, only one copy in the assembly ------
    # refgene6 is in the genome and the annotation. Copy B exists only as transcripts
    # (~12% of codons swapped), standing in for a paralog the assembly collapsed or the
    # annotation missed. Both should map to refgene6 and be reported as one group of two.
    par = place(Gene(rng, "refgene6", n_codons=300), "chr1", "+")
    annotated.append(par)
    gff += gff_lines(par)
    emit_reads(par, 5, "reference")

    cdsA = par.mrna[par.cds_start - 1:par.cds_end]
    cdsB = diverge_cds(rng, cdsA, 0.12)
    mrnaB = (par.mrna[:par.cds_start - 1] + cdsB + par.mrna[par.cds_end:])
    for i in range(5):
        reads.append((f"paralogB_read{i}", mutate(rng, mrnaB, 0.002)))
    ident = 100.0 * sum(a == b for a, b in zip(cdsA, cdsB)) / len(cdsA)
    key.append(("paralogB", f"reference (paralog of refgene6, {ident:.1f}% CDS identity, "
                            f"not in the assembly)"))

    add_spacer("chr1", 1500)
    add_spacer("chr2", 1500)

    with open(out / "genome.fa", "w") as fh:
        for c, parts in chroms.items():
            s = "".join(parts)
            fh.write(f">{c}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    with open(out / "reference.gff3", "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write("\n".join(gff) + "\n")

    rng.shuffle(reads)
    with open(out / "flnc.fa", "w") as fh:
        for n, s in reads:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    with open(out / "flnc.fq", "w") as fh:
        for n, s in reads:
            fh.write(f"@{n}\n{s}\n+\n{'I' * len(s)}\n")

    write_pacbio_bam(reads, out / "flnc.bam")

    with open(out / "expected_classes.tsv", "w") as fh:
        fh.write("locus\texpected_class\n")
        for n, c in key:
            fh.write(f"{n}\t{c}\n")

    print(f"wrote {out}/  genome.fa reference.gff3 flnc.fa flnc.bam ({len(reads)} reads)")


if __name__ == "__main__":
    main()
