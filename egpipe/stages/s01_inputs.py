"""Stage 01 - build the tr2aacds input fasta.

Two evidence tracks are concatenated into one file:
  1. FLNC reads (renamed to short safe IDs; identical sequences optionally collapsed)
  2. reference CDS extracted from the genome + annotation with gffread

Renaming matters: PacBio read names contain '/' and are long, which breaks
several evigene perl scripts and the file names they derive.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from ..utils import LOG, die, iter_fasta, require_tools, run, write_fasta

STAGE = "01_inputs"


def _check_pacbio_bam(bam: Path) -> None:
    """isoseq cluster2 hangs (no error, no output) on a BAM without a PacBio @RG.

    uBAMs made by `samtools import` or similar tools have no read group carrying
    PL:PACBIO / DS:READTYPE=..., so fail fast with a message instead of hanging.
    """
    import pysam
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as fh:
        rgs = fh.header.to_dict().get("RG", [])
    ok = any(rg.get("PL", "").upper() == "PACBIO" and "READTYPE" in rg.get("DS", "")
             for rg in rgs)
    if not ok:
        die(
            f"{bam} has no PacBio read group (PL:PACBIO with DS:READTYPE=...).\n"
            "  isoseq cluster2 hangs silently on such files. Either use the FLNC BAM as\n"
            "  delivered by SMRT Link, or set tracks.cluster.enabled: false to feed the\n"
            "  raw reads to tr2aacds instead."
        )


# Filesystem types that live in RAM. `rootfs` matters as much as `tmpfs` here:
# diskless compute nodes boot with their root in memory, so /tmp on such a node is
# RAM and anything written there is charged to the job's memory cgroup.
RAM_BACKED = ("tmpfs", "ramfs", "rootfs", "initramfs")


def _fs_of(path: Path) -> tuple[str, str]:
    """(filesystem type, mount point) for the closest mount point above path."""
    best = ("", "")
    try:
        real = str(path.resolve())
        with open("/proc/mounts") as fh:
            for line in fh:
                f = line.split()
                if len(f) < 3:
                    continue
                mnt, fstype = f[1], f[2]
                if (real == mnt or real.startswith(mnt.rstrip("/") + "/")) \
                        and len(mnt) >= len(best[1]):
                    best = (fstype, mnt)
    except OSError:
        pass
    return best


def check_scratch(cfg, need_bytes: int = 0) -> None:
    """Fail early if scratch is RAM-backed or too small.

    cluster2 spools a temporary BAM of everything it writes. Sent to a tmpfs /tmp
    those pages are charged to the job's memory cgroup but never show up in MaxRSS,
    so the job is OOM-killed while `sacct` reports a tiny peak RSS - which reads as
    anything but a memory problem.
    """
    import os
    cfg.tmpdir.mkdir(parents=True, exist_ok=True)
    fstype, mnt = _fs_of(cfg.tmpdir)
    st = os.statvfs(cfg.tmpdir)
    free = st.f_bavail * st.f_frsize
    LOG.info("scratch: %s (%s on %s, %.1f GB free)",
             cfg.tmpdir, fstype or "unknown fs", mnt or "?", free / 1e9)
    if fstype in RAM_BACKED:
        die(f"resources.tmpdir ({cfg.tmpdir}) is on {fstype}, which is RAM.\n"
            "  cluster2's temporary BAM would be charged to the job's memory cgroup\n"
            "  and the job will be OOM-killed with a small MaxRSS. Point\n"
            "  resources.tmpdir at a real filesystem.")
    if need_bytes and free < need_bytes:
        LOG.warning("scratch has %.1f GB free but the input BAMs total %.1f GB; "
                    "cluster2 spools a temporary copy", free / 1e9, need_bytes / 1e9)


def _stage_env(cfg) -> dict:
    """Keep scratch on the analysis filesystem, not the compute node.

    cluster2 spools a temporary BAM of every record it is about to write, and it
    picks the directory from TMPDIR (verified: with TMPDIR set, tmp.<pid>.bam.tmp
    is created there rather than in /tmp). On a cluster whose nodes have no local
    scratch, letting it fall back to /tmp gives

        [pbbam] BAM writer ERROR: could not write record
          file: /tmp/tmp.NNNN.bam.tmp
          reason: Resource temporarily unavailable

    part-way through a large library. samtools honours TMPDIR too.
    """
    cfg.tmpdir.mkdir(parents=True, exist_ok=True)
    return {"TMPDIR": str(cfg.tmpdir), "TMP": str(cfg.tmpdir), "TEMP": str(cfg.tmpdir)}


def _cluster2(cfg, work: Path) -> Path:
    """Run isoseq cluster2 on the FLNC BAM(s); return the transcripts BAM."""
    require_tools(["isoseq"])
    ccfg = cfg.d["tracks"]["cluster"]
    bams = cfg.flnc_bam
    if not bams:
        die("tracks.cluster.enabled requires inputs.flnc_bam "
            "(isoseq cluster2 cannot read fasta)")
    for b in bams:
        _check_pacbio_bam(b)

    if len(bams) == 1:
        src = str(bams[0])
    else:
        fofn = work / "flnc.fofn"
        fofn.write_text("\n".join(str(b) for b in bams) + "\n")
        src = str(fofn)

    check_scratch(cfg, sum(b.stat().st_size for b in bams))
    out_bam = work / "clustered.bam"
    threads = int(ccfg.get("threads") or cfg.threads)
    cmd = ["isoseq", "cluster2", "-j", str(threads), "--log-level", "INFO",
           "--log-file", str(work / "cluster2.log")]
    if ccfg.get("singletons", True):
        cmd.append("--singletons")
    # --sort-threads is documented as "open file handles per BAM" and defaults to -j,
    # so it is the one lever cluster2 exposes that trades speed for memory. There is no
    # --max-memory option.
    st = ccfg.get("sort_threads")
    if st:
        cmd += ["--sort-threads", str(int(st))]
    cmd += [str(a) for a in (ccfg.get("extra_args") or [])]
    cmd += [src, str(out_bam)]
    LOG.info("cluster2 scratch (TMPDIR): %s", cfg.tmpdir)
    run(cmd, env=_stage_env(cfg))
    if not out_bam.exists():
        die(f"isoseq cluster2 produced no {out_bam}")
    return out_bam


def _clustered_bam_to_fasta(bam: Path, out_fa: Path, map_tsv: Path, prefix: str) -> int:
    """Write renamed transcripts plus the cluster size (is tag) as dupcount."""
    import pysam
    n = 0
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as fh, \
            open(out_fa, "w") as out, open(map_tsv, "w") as mp:
        mp.write("new_id\toriginal_id\tlength\tdupcount\n")
        for rec in fh:
            seq = rec.query_sequence
            if not seq:
                continue
            n += 1
            new = f"{prefix}_{n:07d}"
            write_fasta(out, new, seq)
            size = rec.get_tag("is") if rec.has_tag("is") else 1
            mp.write(f"{new}\t{rec.query_name}\t{len(seq)}\t{size}\n")
    return n


def _bam_to_fasta(bams, out_fa: Path, threads: int, env: dict | None = None) -> None:
    import os
    full = {**os.environ, **(env or {})}
    with open(out_fa, "w") as fh:
        for bam in bams:
            LOG.info("converting %s -> fasta", bam)
            proc = subprocess.Popen(
                ["samtools", "fasta", "-@", str(max(1, threads - 1)), "-0", "-", str(bam)],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, env=full,
            )
            for line in proc.stdout:
                fh.write(line)
            proc.wait()
            if proc.returncode != 0:
                raise SystemExit(f"samtools fasta failed on {bam}")


def _rename(in_fa: Path, out_fa: Path, map_tsv: Path, prefix: str,
            dedup: bool = True) -> int:
    """Rename to short safe IDs; optionally collapse exactly identical sequences.

    Collapsing is done here rather than with `seqkit rmdup` so the number of
    identical reads behind each retained sequence is preserved (dupcount column)
    - that is the only abundance signal that survives into the evigene run.
    """
    n = 0
    seen: dict[str, str] = {}
    counts: dict[str, int] = {}
    with open(out_fa, "w") as out, open(map_tsv, "w") as mp:
        mp.write("new_id\toriginal_id\tlength\tdupcount\n")
        rows: list[tuple[str, str, int]] = []
        for name, _desc, seq in iter_fasta(in_fa):
            if dedup:
                key = seq.upper()
                if key in seen:
                    counts[seen[key]] += 1
                    continue
                n += 1
                new = f"{prefix}_{n:07d}"
                seen[key] = new
                counts[new] = 1
            else:
                n += 1
                new = f"{prefix}_{n:07d}"
                counts[new] = 1
            write_fasta(out, new, seq)
            rows.append((new, name, len(seq)))
        for new, name, ln in rows:
            mp.write(f"{new}\t{name}\t{ln}\t{counts[new]}\n")
    return n


def _ref_sequences(cfg, work: Path) -> tuple[Path, Path]:
    """Always extract reference CDS + protein: needed for the track AND for the
    post-hoc appendix of reference genes the reads never recovered."""
    cds = work / "ref_cds.raw.fa"
    prot = work / "ref_prot.raw.fa"
    if not cds.exists() or not prot.exists():
        run(["gffread", "-x", str(cds), "-y", str(prot),
             "-g", str(cfg.genome), str(cfg.annotation)], env=_stage_env(cfg))
    return cds, prot


def _ref_cds_track(cfg, work: Path) -> tuple[Path, Path, int]:
    raw, _prot = _ref_sequences(cfg, work)
    prefix = cfg.d["tracks"]["ref_cds_prefix"]
    out_fa = work / "ref_cds.fa"
    map_tsv = work / "ref_cds.idmap.tsv"
    n = 0
    with open(out_fa, "w") as out, open(map_tsv, "w") as mp:
        mp.write("new_id\tref_transcript_id\tlength\n")
        for name, _desc, seq in iter_fasta(raw):
            if not seq:
                continue
            n += 1
            new = f"{prefix}_{n:06d}"
            write_fasta(out, new, seq)
            mp.write(f"{new}\t{name}\t{len(seq)}\n")
    return out_fa, map_tsv, n


def run_stage(cfg) -> Path:
    require_tools(["samtools", "gffread"])
    work = cfg.dir("01_inputs")
    prefix = cfg.d["tracks"]["flnc_prefix"]

    # --- FLNC track ----------------------------------------------------------
    cluster_cfg = cfg.d["tracks"].get("cluster", {})
    use_cluster = bool(cluster_cfg.get("enabled", False))
    parts: list[Path] = []
    n_flnc = n_raw = 0

    if use_cluster:
        clustered = _cluster2(cfg, work)
        flnc_fa = work / "flnc.fa"
        n_flnc = _clustered_bam_to_fasta(clustered, flnc_fa, work / "flnc.idmap.tsv", prefix)
        LOG.info("isoseq cluster2 track: %d transcripts -> %s", n_flnc, flnc_fa)
        parts.append(flnc_fa)

        if cfg.d["tracks"].get("keep_raw_flnc_track", False):
            raw = work / "flnc.raw.fa"
            _bam_to_fasta(cfg.flnc_bam, raw, cfg.threads, _stage_env(cfg))
            raw_fa = work / "flnc_raw.fa"
            n_raw = _rename(raw, raw_fa, work / "flnc_raw.idmap.tsv",
                            prefix + "raw", dedup=True)
            LOG.info("raw FLNC track (extra): %d sequences", n_raw)
            parts.append(raw_fa)
    else:
        check_scratch(cfg)
        flnc_raw = work / "flnc.raw.fa"
        if cfg.flnc_bam:
            _bam_to_fasta(cfg.flnc_bam, flnc_raw, cfg.threads, _stage_env(cfg))
        if cfg.flnc_fasta:
            mode = "a" if cfg.flnc_bam else "w"
            with open(flnc_raw, mode) as out:
                for fa in cfg.flnc_fasta:
                    for name, _d, seq in iter_fasta(fa):
                        write_fasta(out, name, seq)
        flnc_fa = work / "flnc.fa"
        n_flnc = _rename(flnc_raw, flnc_fa, work / "flnc.idmap.tsv", prefix,
                         dedup=bool(cfg.d["tracks"]["dedup_identical"]))
        LOG.info("raw FLNC track: %d sequences -> %s", n_flnc, flnc_fa)
        parts.append(flnc_fa)

    # --- reference sequences -------------------------------------------------
    # extracted unconditionally; only added as a tr2aacds evidence track on request
    _ref_sequences(cfg, work)
    n_ref = 0
    if cfg.d["tracks"]["include_reference_cds"]:
        ref_fa, _map, n_ref = _ref_cds_track(cfg, work)
        parts.append(ref_fa)
        LOG.info("reference CDS track: %d sequences -> %s", n_ref, ref_fa)
    else:
        LOG.info("reference CDS not added as a tr2aacds track "
                 "(unrecovered reference genes are appended after classification)")

    merged = work / f"{cfg.evg_name}.fa"
    with open(merged, "w") as out:
        for p in parts:
            with open(p) as fh:
                for line in fh:
                    out.write(line)

    label = "isoseq_cluster2" if use_cluster else "flnc_reads"
    total = n_flnc + n_raw + n_ref
    (work / "input_counts.tsv").write_text(
        "track\tn_sequences\n"
        f"{label}\t{n_flnc}\nflnc_raw_extra\t{n_raw}\n"
        f"reference_cds\t{n_ref}\ntotal\t{total}\n"
    )
    LOG.info("tr2aacds input: %s (%d sequences)", merged, total)
    return merged
