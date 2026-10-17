"""Reference annotation model: CDS exons and gene spans, indexed per sequence.

Only coding features are used - the pipeline is CDS-focused by design.
Supports GFF3 (ID/Parent) and GTF (gene_id/transcript_id).
"""
from __future__ import annotations

import gzip
import re
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

from .utils import LOG, die

ATTR_GFF3 = re.compile(r"([^=;\s]+)=([^;]*)")
ATTR_GTF = re.compile(r'(\S+)\s+"([^"]*)"')


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _sniff(path: Path, declared: str) -> str:
    if declared in ("gff3", "gtf"):
        return declared
    if str(path).lower().endswith((".gtf", ".gtf.gz")):
        return "gtf"
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            attrs = line.rstrip("\n").split("\t")[-1]
            if "gene_id " in attrs and "=" not in attrs.split(";")[0]:
                return "gtf"
            return "gff3"
    return "gff3"


class RefAnnot:
    """CDS exons and gene spans, with sorted per-sequence indexes for overlap queries."""

    def __init__(self):
        self.cds: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)   # chrom -> (s,e,strand,gene)
        self.genes: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)  # chrom -> (s,e,strand,gene)
        self.tx2gene: dict[str, str] = {}
        self._cds_starts: dict[str, list[int]] = {}
        self._gene_starts: dict[str, list[int]] = {}
        self._cds_maxend: dict[str, list[int]] = {}
        self._gene_maxend: dict[str, list[int]] = {}

    # -------------------------------------------------------------- building
    def finalize(self) -> None:
        for store, starts, maxend in (
            (self.cds, self._cds_starts, self._cds_maxend),
            (self.genes, self._gene_starts, self._gene_maxend),
        ):
            for chrom, items in store.items():
                items.sort(key=lambda t: (t[0], t[1]))
                starts[chrom] = [t[0] for t in items]
                run = []
                cur = 0
                for t in items:
                    cur = max(cur, t[1])
                    run.append(cur)
                maxend[chrom] = run

    def _query(self, store, starts, maxend, chrom, s, e):
        items = store.get(chrom)
        if not items:
            return []
        st = starts[chrom]
        me = maxend[chrom]
        # candidates: all intervals whose start <= e; prune with running max end
        hi = bisect_left(st, e + 1)
        res = []
        i = hi - 1
        while i >= 0:
            if me[i] < s:
                break
            a, b, strand, gid = items[i]
            if a <= e and b >= s:
                res.append(items[i])
            i -= 1
        return res

    def cds_overlaps(self, chrom, s, e):
        return self._query(self.cds, self._cds_starts, self._cds_maxend, chrom, s, e)

    def gene_overlaps(self, chrom, s, e):
        return self._query(self.genes, self._gene_starts, self._gene_maxend, chrom, s, e)

    @property
    def n_cds(self) -> int:
        return sum(len(v) for v in self.cds.values())

    @property
    def n_genes(self) -> int:
        return sum(len(v) for v in self.genes.values())


def load(path: Path, fmt: str = "auto") -> RefAnnot:
    path = Path(path)
    fmt = _sniff(path, fmt)
    LOG.info("parsing reference annotation (%s): %s", fmt, path)
    ra = RefAnnot()

    tx_parent: dict[str, str] = {}      # transcript ID -> gene ID
    gene_span: dict[str, list] = {}     # gene ID -> [chrom, start, end, strand]
    cds_rows: list[tuple[str, int, int, str, str]] = []   # chrom,s,e,strand,parent_tx

    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _src, feat, s, e, _sc, strand, _ph, attrs = f[:9]
            s, e = int(s), int(e)
            featl = feat.lower()
            if fmt == "gff3":
                kv = dict(ATTR_GFF3.findall(attrs))
                fid = kv.get("ID")
                parent = (kv.get("Parent") or "").split(",")[0]
            else:
                kv = dict(ATTR_GTF.findall(attrs))
                fid = kv.get("transcript_id")
                parent = kv.get("gene_id")

            if featl == "gene":
                gid = fid or parent
                if gid:
                    gene_span[gid] = [chrom, s, e, strand]
            elif featl in ("mrna", "transcript") or featl.endswith("_gene_segment"):
                gid = parent or fid
                if fid:
                    tx_parent[fid] = gid or fid
                if gid and gid not in gene_span:
                    gene_span[gid] = [chrom, s, e, strand]
            elif featl == "cds":
                tx = parent or fid
                if fmt == "gtf":
                    tx = kv.get("transcript_id") or tx
                    if tx and kv.get("gene_id"):
                        tx_parent.setdefault(tx, kv["gene_id"])
                cds_rows.append((chrom, s, e, strand, tx or "NA"))

    if not cds_rows:
        die(f"no CDS features found in {path} - the pipeline is CDS-only and needs them")

    coding_genes: set[str] = set()
    for chrom, s, e, strand, tx in cds_rows:
        gid = tx_parent.get(tx, tx)
        coding_genes.add(gid)
        ra.cds[chrom].append((s, e, strand, gid))
        ra.tx2gene[tx] = gid
        span = gene_span.get(gid)
        if span is None:
            gene_span[gid] = [chrom, s, e, strand]
        else:
            span[1] = min(span[1], s)
            span[2] = max(span[2], e)

    # gene spans are restricted to CODING genes: this pipeline is CDS-only, so a
    # transcript falling inside a non-coding locus still counts as intergenic.
    for gid, (chrom, s, e, strand) in gene_span.items():
        if gid in coding_genes:
            ra.genes[chrom].append((s, e, strand, gid))

    ra.finalize()
    LOG.info("reference: %d genes, %d CDS exons across %d sequences",
             ra.n_genes, ra.n_cds, len(ra.cds))
    return ra
