"""Stage 03 - harvest the evigene okayset into a clean transcriptome.

Produces the final transcript / CDS / protein fasta plus an ORF table holding
the CDS offsets *in transcript coordinates* (this is what the GFF is built on).

evigene fasta headers look like:
  >Em2EVm000012t1 type=mRNA; aalen=616,59%,complete; clen=3115; strand=+;
   offs=455-2305; oid=flnc_0001234,refcds_000512; ...
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

from ..utils import LOG, iter_fasta, revcomp, write_fasta

STAGE = "03_harvest"

KV = re.compile(r"(\w+)=([^;]+)")


def parse_header(desc: str) -> dict:
    d = {}
    for k, v in KV.findall(desc or ""):
        d[k] = v.strip()
    return d


def parse_orf(desc: str, seqlen: int) -> dict:
    """Return cds_start/cds_end (1-based, always start<=end), cds_strand, aalen, completeness."""
    kv = parse_header(desc)
    out = {
        "clen": int(kv.get("clen", seqlen) or seqlen),
        "cds_start": None, "cds_end": None, "cds_strand": "+",
        "aa_len": None, "pct_cds": None, "completeness": None,
        "oids": kv.get("oid", ""), "evg_type": kv.get("type", ""),
    }
    aalen = kv.get("aalen")
    if aalen:
        bits = aalen.split(",")
        try:
            out["aa_len"] = int(bits[0])
        except ValueError:
            pass
        if len(bits) > 1:
            out["pct_cds"] = bits[1].rstrip("%")
        if len(bits) > 2:
            out["completeness"] = bits[2]
    offs = kv.get("offs")
    if offs:
        m = re.match(r"^(-?\d+)-(-?\d+)", offs)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            strand = "+"
            if a > b:
                a, b = b, a
                strand = "-"
            if kv.get("strand") == "-" and strand == "+" and kv.get("revc") != "1":
                strand = "-"
            out["cds_start"], out["cds_end"], out["cds_strand"] = a, b, strand
    return out


def _find_seqsets(evgdir: Path, include_alt: bool,
                  include_cull: bool = False) -> dict[str, list[Path]]:
    """Locate mrna/cds/aa fasta in okayset, handling both evigene layouts.

    evigene may also write a `<name>.cull.*` set alongside `<name>.okay.*`. Culls
    sit between okay and drop - trclass2pubset filters them with `/(drop|cull)/` -
    so they are excluded by default and included only on request.
    """
    okay = evgdir / "okayset"
    sets: dict[str, list[Path]] = {}
    for ext in ("mrna", "cds", "aa"):
        pub = sorted(
            p for p in okay.glob(f"*.{ext}")
            if not re.search(r"\.(okay|okalt|okreor|cull)\.", p.name)
        )
        pub += sorted(okay.glob(f"*.{ext}_pub.fa"))
        if pub:
            sets[ext] = pub[:1]
            continue
        staged = sorted(okay.glob(f"*.okay.{ext}"))
        if include_alt:
            staged += sorted(okay.glob(f"*.okalt.{ext}"))
        if include_cull:
            staged += sorted(okay.glob(f"*.cull.{ext}"))
        if not staged:
            alt_dir = evgdir / "okayset1st"
            staged = sorted(alt_dir.glob(f"*.okay.{ext}"))
            if include_alt:
                staged += sorted(alt_dir.glob(f"*.okalt.{ext}"))
        sets[ext] = staged
    return sets


def read_pubids(evgdir: Path) -> dict[str, dict]:
    """pubid -> {orig_id, gene, altnum, class, aaqual, notes}."""
    out: dict[str, dict] = {}
    cands = list((evgdir / "okayset").glob("*.pubids")) + list((evgdir / "publicset").glob("*.pubids"))
    cands = [c for c in cands if not c.name.endswith(".old")]
    if not cands:
        return out
    with open(cands[0]) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            out[f[0]] = {
                "orig_id": f[1], "locus": f[2], "alt_num": f[3], "evg_class": f[4],
                "aa_qual": f[5] if len(f) > 5 else "",
                "notes": f[7] if len(f) > 7 else "",
                "oids": f[8] if len(f) > 8 else "",
            }
    LOG.info("pubids: %d entries from %s", len(out), cands[0].name)
    return out


def read_trclass(evgdir: Path) -> dict[str, dict]:
    """orig_id -> {keep, evg_class, match_id, pct}."""
    out: dict[str, dict] = {}
    cands = sorted(evgdir.glob("*.trclass")) + sorted((evgdir / "okayset").glob("*.trclass"))
    if not cands:
        return out
    with open(cands[0]) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            out[f[0]] = {"keep": f[1], "evg_class": f[2], "match_id": f[3],
                         "pct": f[4] if len(f) > 4 else ""}
    LOG.info("trclass: %d entries from %s", len(out), cands[0].name)
    return out


def load_oid_chain(cfg) -> dict[str, str]:
    """round-1 public ID -> original input IDs, written by the chunked stage 02."""
    p = cfg.outdir / "02_evigene" / "oid_chain.tsv"
    out: dict[str, str] = {}
    if not p.exists():
        return out
    with open(p) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = f[1]
    LOG.info("resolving oids through the round-1 chain (%d entries)", len(out))
    return out


def _resolve_oids(oids: str, chain: dict[str, str]) -> str:
    if not chain or not oids:
        return oids
    out: list[str] = []
    for o in str(oids).split(","):
        o = o.strip()
        if not o:
            continue
        out.extend(x for x in chain.get(o, o).split(",") if x)
    seen, uniq = set(), []
    for o in out:
        if o not in seen:
            seen.add(o)
            uniq.append(o)
    return ",".join(uniq)


def run_stage(cfg, evgdir: Path) -> dict:
    include_alt = bool(cfg.d["transcriptome"]["include_alternates"])
    min_cds = int(cfg.d["transcriptome"]["min_cds_len"])
    out = cfg.dir("03_transcriptome")

    include_cull = bool(cfg.d["transcriptome"].get("include_cull", False))
    sets = _find_seqsets(evgdir, include_alt, include_cull)
    if not sets.get("mrna"):
        raise SystemExit(f"no okayset mRNA fasta found under {evgdir}/okayset")
    LOG.info("okayset mRNA files: %s", ", ".join(p.name for p in sets["mrna"]))
    if not include_cull:
        culled = sorted((evgdir / "okayset").glob("*.cull.mrna"))
        n_cull = sum(sum(1 for line in open(c) if line.startswith(">")) for c in culled)
        if n_cull:
            LOG.info("excluding %d transcripts in %s (evigene 'cull' class - it groups "
                     "these with drops; set transcriptome.include_cull: true to keep them)",
                     n_cull, ", ".join(c.name for c in culled))

    pubids = read_pubids(evgdir)
    trclass = read_trclass(evgdir)
    chain = load_oid_chain(cfg)

    tx_fa = out / f"{cfg.sample}.transcripts.fa"
    orf_tsv = out / f"{cfg.sample}.orf_table.tsv"
    ref_prefix = cfg.d["tracks"]["ref_cds_prefix"]

    # ---- pass 1: read headers, build records, decide what to keep ----------
    cand: dict[str, dict] = {}
    n_in = n_noorf = n_short = 0
    for path in sets["mrna"]:
        for name, desc, seq in iter_fasta(path):
            n_in += 1
            info = parse_orf(desc, len(seq))
            if info["cds_start"] is None:
                n_noorf += 1
                continue
            cds_len = info["cds_end"] - info["cds_start"] + 1
            if min_cds and cds_len < min_cds:
                n_short += 1
                continue
            if name in cand:
                continue
            pub = pubids.get(name, {})
            oids = _resolve_oids(info["oids"] or pub.get("oids", ""), chain)
            cand[name] = {
                "transcript_id": name, "length": len(seq),
                "cds_start": info["cds_start"], "cds_end": info["cds_end"],
                "cds_strand": info["cds_strand"], "aa_len": info["aa_len"],
                "pct_cds": info["pct_cds"], "completeness": info["completeness"],
                "locus": pub.get("locus", ""), "alt_num": pub.get("alt_num", ""),
                "evg_class": pub.get("evg_class", "") or info["evg_type"],
                "oids": oids,
                "_cds_len": cds_len,
                "_refcds_only": bool(oids) and all(
                    o.strip().startswith(ref_prefix) for o in str(oids).split(",") if o.strip()),
            }

    n_alt_skipped = 0
    if include_alt:
        keep_ids = set(cand)
    else:
        # one transcript per evigene locus. Prefer an FLNC-derived sequence over a
        # reference-CDS placeholder even when evigene made the placeholder the main.
        by_locus: dict[str, list[dict]] = {}
        for rec in cand.values():
            by_locus.setdefault(rec["locus"] or rec["transcript_id"], []).append(rec)
        keep_ids = set()
        for _loc, recs in by_locus.items():
            flnc = [r for r in recs if not r["_refcds_only"]]
            pool = flnc or recs
            main = [r for r in pool if not str(r["evg_class"]).startswith("alt")]
            best = max(main or pool, key=lambda r: (r["_cds_len"], r["length"]))
            keep_ids.add(best["transcript_id"])
        n_alt_skipped = len(cand) - len(keep_ids)

    # ---- pass 2: write the kept sequences and the ORF table ----------------
    # Every transcript is stored ORF-sense. evigene normally does this itself
    # (revc=1 in the header), but it can emit a reversed ORF offset without having
    # flipped the sequence. Fixing it here means the transcriptome is uniformly
    # sense-oriented, which is what makes `minimap2 -u f` valid at stage 04.
    kept: dict[str, dict] = {}
    n_flipped = 0
    with open(tx_fa, "w") as fh, open(orf_tsv, "w", newline="") as th:
        w = csv.writer(th, delimiter="\t")
        w.writerow(["transcript_id", "length", "cds_start", "cds_end", "cds_strand",
                    "aa_len", "pct_cds", "completeness", "locus", "alt_num",
                    "evg_class", "oids"])
        for path in sets["mrna"]:
            for name, desc, seq in iter_fasta(path):
                if name not in keep_ids or name in kept:
                    continue
                rec = {k: v for k, v in cand[name].items() if not k.startswith("_")}
                if rec["cds_strand"] == "-":
                    L = len(seq)
                    seq = revcomp(seq)
                    rec["cds_start"], rec["cds_end"] = (
                        L - rec["cds_end"] + 1, L - rec["cds_start"] + 1)
                    rec["cds_strand"] = "+"
                    desc = re.sub(r"offs=[^;]*;?",
                                  f"offs={rec['cds_start']}-{rec['cds_end']}; reoriented=1;",
                                  desc or "")
                    n_flipped += 1
                write_fasta(fh, name, seq, desc)
                kept[name] = rec
                w.writerow([rec[k] for k in
                            ("transcript_id", "length", "cds_start", "cds_end", "cds_strand",
                             "aa_len", "pct_cds", "completeness", "locus", "alt_num",
                             "evg_class", "oids")])

    # subset cds / aa fasta to the kept ids
    for ext, label in (("cds", "cds"), ("aa", "protein")):
        if not sets.get(ext):
            continue
        dest = out / f"{cfg.sample}.{label}.fa"
        n = 0
        with open(dest, "w") as fh:
            for path in sets[ext]:
                for name, desc, seq in iter_fasta(path):
                    base = name[:-6] if name.endswith("utrorf") else name
                    if name in kept or base in kept:
                        write_fasta(fh, name, seq, desc)
                        n += 1
        LOG.info("wrote %s (%d seqs)", dest.name, n)

    LOG.info("harvest: %d input, %d kept, %d no-ORF, %d below min_cds_len, %d alternates skipped",
             n_in, len(kept), n_noorf, n_short, n_alt_skipped)
    if n_flipped:
        LOG.info("reoriented %d transcripts to ORF sense (evigene left them reversed); "
                 "the transcriptome is now uniformly sense-oriented", n_flipped)
    return {"transcripts_fa": tx_fa, "orf_table": orf_tsv, "records": kept,
            "pubids": pubids, "trclass": trclass}
