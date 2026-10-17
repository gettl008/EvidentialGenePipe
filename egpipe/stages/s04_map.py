"""Stage 04 - splice-aware alignment of the final transcripts to the genome."""
from __future__ import annotations

from pathlib import Path

from ..utils import LOG, require_tools, run

STAGE = "04_map"


def run_stage(cfg, transcripts_fa: Path) -> Path:
    require_tools(["minimap2", "samtools"])
    out = cfg.dir("04_mapping")
    bam = out / f"{cfg.sample}.transcripts.genome.bam"

    mp = cfg.d["mapping"]
    cmd = [
        "minimap2", "-a", "-x", str(mp["preset"]), "-u", str(mp["strand_mode"]),
        "-t", str(cfg.threads), "--cs=long", "-C5",
    ]
    if not mp.get("secondary", False):
        cmd += ["--secondary=no"]
    cmd += [str(a) for a in (mp.get("extra_args") or [])]
    cmd += [str(cfg.genome), str(transcripts_fa)]

    sam = out / f"{cfg.sample}.transcripts.genome.sam"
    env = {"TMPDIR": str(cfg.tmpdir), "TMP": str(cfg.tmpdir), "TEMP": str(cfg.tmpdir)}
    cfg.tmpdir.mkdir(parents=True, exist_ok=True)
    run(cmd, stdout=sam, env=env)
    run(["samtools", "sort", "-@", str(max(1, cfg.threads - 1)),
         "-T", str(cfg.tmpdir / "mmsort"), "-o", str(bam), str(sam)], env=env)
    run(["samtools", "index", "-@", str(max(1, cfg.threads - 1)), str(bam)])
    sam.unlink(missing_ok=True)

    stats = out / f"{cfg.sample}.flagstat.txt"
    run(["samtools", "flagstat", str(bam)], stdout=stats)
    LOG.info("alignment: %s", bam)
    return bam
