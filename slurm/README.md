# SLURM usage on RDAC

Generate per-stage scripts (nothing is hardcoded; they are rendered from the config):

    python -m egpipe sbatch --config config/Em2.yaml --partition defq --conda-env evigene_pipe

That writes `<outdir>/slurm/Em2.stage01.sbatch` ... `stage06.sbatch` plus
`Em2.submit_all.sh`, which submits all six with `afterok` dependencies.

Resource notes for RDAC:

* `defq` nodes are 32 CPU / ~83 GB. The stage-02 script uses `resources.threads`
  and `resources.max_mem_mb` from the config, so keep `max_mem_mb` at or below
  ~76000 for defq.
* Stage 02 (tr2aacds) is the long one; its slowest step is the all-vs-all
  self-blastn of CDS. With raw FLNC input this is the stage to watch. If it
  needs more than 83 GB, re-render with `--partition bigmem` and raise
  `resources.max_mem_mb`.
* `TmpDisk=0` on all RDAC nodes, so `resources.tmpdir` should point at a
  filesystem path (default `<outdir>/tmp`), never `/tmp` scratch.
