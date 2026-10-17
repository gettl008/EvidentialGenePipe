"""Shared helpers: logging, subprocess, checkpoints, fasta IO."""
from __future__ import annotations

import gzip
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterator, Sequence

LOG = logging.getLogger("egpipe")


def setup_logging(logfile: Path | None = None, verbose: bool = True) -> None:
    LOG.handlers.clear()
    LOG.setLevel(logging.DEBUG if verbose else logging.INFO)
    fmt = logging.Formatter("[%(asctime)s] %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    LOG.addHandler(sh)
    if logfile is not None:
        logfile.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(logfile)
        fh.setFormatter(fmt)
        LOG.addHandler(fh)


def die(msg: str) -> None:
    LOG.error(msg)
    raise SystemExit(1)


def require_tools(tools: Sequence[str]) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        die(
            "missing executables on PATH: %s\n"
            "  activate the pipeline env:  conda activate evigene_pipe" % ", ".join(missing)
        )


def run(cmd: Sequence[str] | str, cwd: Path | None = None, env: dict | None = None,
        stdout: Path | None = None, check: bool = True) -> int:
    """Run a command, streaming stderr to the log."""
    shell = isinstance(cmd, str)
    pretty = cmd if shell else " ".join(str(c) for c in cmd)
    LOG.info("RUN  %s%s", pretty, f"  > {stdout}" if stdout else "")
    t0 = time.time()
    full_env = dict(os.environ)
    if env:
        full_env.update({k: str(v) for k, v in env.items()})
    out_fh = open(stdout, "w") if stdout else None
    try:
        proc = subprocess.run(
            cmd if shell else [str(c) for c in cmd],
            shell=shell, cwd=str(cwd) if cwd else None, env=full_env,
            stdout=out_fh if out_fh else None,
        )
    finally:
        if out_fh:
            out_fh.close()
    dt = time.time() - t0
    LOG.info("DONE rc=%d in %.1fs", proc.returncode, dt)
    if check and proc.returncode != 0:
        die(f"command failed (rc={proc.returncode}): {pretty}")
    return proc.returncode


# ---------------------------------------------------------------- checkpoints
def stamp_path(outdir: Path, stage: str) -> Path:
    return outdir / "checkpoints" / f"{stage}.done"


def is_done(outdir: Path, stage: str) -> bool:
    return stamp_path(outdir, stage).exists()


def mark_done(outdir: Path, stage: str) -> None:
    p = stamp_path(outdir, stage)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(time.strftime("%Y-%m-%d %H:%M:%S\n"))


def clear_done(outdir: Path, stage: str) -> None:
    stamp_path(outdir, stage).unlink(missing_ok=True)


# ---------------------------------------------------------------- fasta
def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def iter_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    """Yield (id, description_rest, sequence)."""
    name = desc = None
    chunks: list[str] = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(chunks)
                head = line[1:].rstrip("\n")
                parts = head.split(None, 1)
                name = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, desc, "".join(chunks)


def write_fasta(fh, name: str, seq: str, desc: str = "", width: int = 60) -> None:
    fh.write(f">{name}{(' ' + desc) if desc else ''}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


REVCOMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(seq: str) -> str:
    return seq.translate(REVCOMP)[::-1]


def first_existing(paths) -> Path | None:
    for p in paths:
        p = Path(p)
        if p.exists() and p.stat().st_size > 0:
            return p
    return None
