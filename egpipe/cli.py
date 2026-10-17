"""Command line entry point.

  python -m egpipe run     --config config/Em2.yaml
  python -m egpipe run     --config config/Em2.yaml --from 04 --force
  python -m egpipe check   --config config/Em2.yaml
  python -m egpipe sbatch  --config config/Em2.yaml   # print/emit SLURM scripts
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .config import load_config
from .stages import s01_inputs, s02_tr2aacds, s03_harvest, s04_map, s05_classify, s06_annotate
from .utils import LOG, clear_done, is_done, mark_done, setup_logging

STAGES = ["01", "02", "03", "04", "05", "06"]
STAGE_NAMES = {
    "01": "inputs (FLNC -> fasta, reference CDS track)",
    "02": "tr2aacds (EvidentialGene)",
    "03": "harvest okayset -> transcriptome + ORF table",
    "04": "minimap2 transcripts -> genome",
    "05": "classify reference vs novel",
    "06": "write transcript-coordinate GFF + deliverables",
}


def _stage_range(start: str, only: str | None) -> list[str]:
    if only:
        return [only]
    i = STAGES.index(start)
    return STAGES[i:]


def cmd_run(args) -> int:
    cfg = load_config(args.config)
    cfg.force = bool(args.force)
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    setup_logging(cfg.outdir / f"{cfg.sample}.egpipe.log")
    LOG.info("sample=%s outdir=%s threads=%d", cfg.sample, cfg.outdir, cfg.threads)

    todo = _stage_range(args.start, args.only)
    LOG.info("stages to run: %s", ", ".join(todo))

    # Outside SLURM nothing caps resources.threads, so a stage launched on a login
    # node runs at full width there. Stages 01, 02 and 04 are the long, wide ones.
    import os
    heavy = [s for s in todo if s in ("01", "02", "04")]
    if heavy and not os.environ.get("SLURM_JOB_ID"):
        LOG.warning("stage(s) %s outside SLURM: %d threads uncapped on this host. "
                    "If this is a login node, stop and submit instead:\n"
                    "    python -m egpipe sbatch --config <config> --partition defq",
                    "/".join(heavy), cfg.threads)

    input_fa = cfg.outdir / "01_inputs" / f"{cfg.evg_name}.fa"
    evgdir = cfg.outdir / "02_evigene"
    for sub in ("final", "round2"):      # chunked run: results live in a subdir
        if (evgdir / sub / "okayset").is_dir():
            evgdir = evgdir / sub
            break
    bam = cfg.outdir / "04_mapping" / f"{cfg.sample}.transcripts.genome.bam"
    harvest = classified = None

    if "01" in todo:
        if is_done(cfg.outdir, "01") and not cfg.force:
            LOG.info("stage 01 already done - skipping (use --force to redo)")
        else:
            clear_done(cfg.outdir, "01")
            input_fa = s01_inputs.run_stage(cfg)
            mark_done(cfg.outdir, "01")

    if "02" in todo:
        if is_done(cfg.outdir, "02") and not cfg.force:
            LOG.info("stage 02 already done - skipping (use --force to redo)")
        else:
            if not input_fa.exists():
                LOG.error("missing %s - run stage 01 first", input_fa)
                return 1
            clear_done(cfg.outdir, "02")
            evgdir = s02_tr2aacds.run_stage(cfg, input_fa)
            mark_done(cfg.outdir, "02")

    # stages 03/05/06 are cheap and always recomputed so the tables stay in sync
    if any(s in todo for s in ("03", "05", "06")):
        harvest = s03_harvest.run_stage(cfg, evgdir)

    # stage 04 needs only the transcript fasta, so `--only 04` can reuse the file
    # stage 03 already wrote instead of recomputing the harvest
    tx_fa = (harvest["transcripts_fa"] if harvest is not None
             else cfg.outdir / "03_transcriptome" / f"{cfg.sample}.transcripts.fa")

    if "04" in todo:
        if not tx_fa.exists():
            LOG.error("missing %s - run stage 03 first (--from 03, or --only 03)", tx_fa)
            return 1
        stale = bam.exists() and tx_fa.stat().st_mtime > bam.stat().st_mtime
        if stale:
            LOG.info("transcript fasta is newer than the alignment - redoing stage 04")
        if is_done(cfg.outdir, "04") and not cfg.force and bam.exists() and not stale:
            LOG.info("stage 04 already done - skipping (use --force to redo)")
        else:
            clear_done(cfg.outdir, "04")
            bam = s04_map.run_stage(cfg, tx_fa)
            mark_done(cfg.outdir, "04")

    if "05" in todo or "06" in todo:
        if harvest is None:
            LOG.error("stage %s needs the harvest from stage 03", "/".join(todo))
            return 1
        if not bam.exists():
            LOG.error("missing %s - run stage 04 first", bam)
            return 1
        classified = s05_classify.run_stage(cfg, harvest, bam)
        mark_done(cfg.outdir, "05")

    if "06" in todo:
        if classified is None:
            LOG.error("stage 06 needs the classification from stage 05")
            return 1
        s06_annotate.run_stage(cfg, harvest, classified, bam)
        mark_done(cfg.outdir, "06")

    LOG.info("finished")
    return 0


def cmd_check(args) -> int:
    setup_logging(None)
    cfg = load_config(args.config)
    from .utils import require_tools
    require_tools(["perl", "fastanrdb", "cd-hit", "cd-hit-est", "blastn", "makeblastdb",
                   "minimap2", "samtools", "gffread"])
    s02_tr2aacds.check_evigene_perl(cfg)
    LOG.info("config OK, all executables present")
    LOG.info("evigene script: %s", cfg.evigene_script)
    for s in STAGES:
        LOG.info("  stage %s  %-55s %s", s, STAGE_NAMES[s],
                 "DONE" if is_done(cfg.outdir, s) else "-")
    print("\ntr2aacds command that will be run:")
    print("  " + " ".join(s02_tr2aacds.build_command(
        cfg, cfg.outdir / "01_inputs" / f"{cfg.evg_name}.fa")))
    return 0


def cmd_sbatch(args) -> int:
    setup_logging(None)
    cfg = load_config(args.config)
    tpl = Path(__file__).resolve().parent.parent / "slurm" / "stage.sbatch.template"
    if not tpl.exists():
        LOG.error("template not found: %s", tpl)
        return 1
    outdir = cfg.dir("slurm")
    repo = Path(__file__).resolve().parent.parent
    written = []
    profiles = {
        # stage 01 runs isoseq cluster2 over every FLNC read, which is the second most
        # memory-hungry step after tr2aacds - the old hardcoded 8 CPU / 32G was wishful
        "01": dict(cpus=cfg.threads, mem=f"{cfg.max_mem_mb}M", time="12:00:00",
                   part=args.partition),
        "02": dict(cpus=cfg.threads, mem=f"{cfg.max_mem_mb}M", time="7-00:00:00", part=args.partition),
        "03": dict(cpus=2, mem="16G", time="04:00:00", part=args.partition),
        "04": dict(cpus=cfg.threads, mem="64G", time="12:00:00", part=args.partition),
        "05": dict(cpus=cfg.threads, mem="48G", time="08:00:00", part=args.partition),
        "06": dict(cpus=2, mem="16G", time="04:00:00", part=args.partition),
    }
    body = tpl.read_text()
    for stage, p in profiles.items():
        txt = (body
               .replace("@JOBNAME@", f"eg{stage}_{cfg.sample}")
               .replace("@PARTITION@", p["part"])
               .replace("@CPUS@", str(p["cpus"]))
               .replace("@MEM@", p["mem"])
               .replace("@TIME@", p["time"])
               .replace("@LOGDIR@", str(outdir))
               .replace("@REPO@", str(repo))
               .replace("@CONFIG@", str(Path(args.config).resolve()))
               .replace("@CONDA_ENV@", args.conda_env)
               .replace("@STAGE@", stage))
        dest = outdir / f"{cfg.sample}.stage{stage}.sbatch"
        dest.write_text(txt)
        written.append(dest)
    chain = outdir / f"{cfg.sample}.submit_all.sh"
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", "cd \"$(dirname \"$0\")\"", "prev=\"\""]
    for stage in STAGES:
        f = f"{cfg.sample}.stage{stage}.sbatch"
        lines.append(f'jid=$(sbatch ${{prev:+--dependency=afterok:$prev}} --parsable {f})')
        lines.append(f'echo "stage {stage} -> $jid"')
        lines.append('prev=$jid')
    chain.write_text("\n".join(lines) + "\n")
    chain.chmod(0o755)
    for w in written:
        LOG.info("wrote %s", w)
    LOG.info("wrote %s  (submits the whole chain with afterok dependencies)", chain)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="egpipe", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run the pipeline")
    r.add_argument("--config", required=True)
    r.add_argument("--from", dest="start", default="01", choices=STAGES,
                   help="start at this stage (default 01)")
    r.add_argument("--only", default=None, choices=STAGES, help="run just this stage")
    r.add_argument("--force", action="store_true", help="ignore checkpoints and redo")
    r.set_defaults(func=cmd_run)

    c = sub.add_parser("check", help="validate config, tools and stage status")
    c.add_argument("--config", required=True)
    c.set_defaults(func=cmd_check)

    s = sub.add_parser("sbatch", help="write SLURM scripts for each stage")
    s.add_argument("--config", required=True)
    s.add_argument("--partition", default="defq")
    s.add_argument("--conda-env", default="evigene_pipe")
    s.set_defaults(func=cmd_sbatch)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
