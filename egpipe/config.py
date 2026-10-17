"""Config loading and validation. Everything the pipeline needs comes from here."""
from __future__ import annotations

import copy
import os
import re
from pathlib import Path
from typing import Any

import yaml

from .utils import LOG, die

DEFAULTS: dict[str, Any] = {
    "sample": None,
    "inputs": {
        "flnc_bam": [],
        "flnc_fasta": [],
        "genome": None,
        "annotation": None,
        "annotation_format": "auto",
    },
    "outdir": None,
    "resources": {"threads": 8, "max_mem_mb": 16000, "tmpdir": None},
    "evigene": {
        "home": None,
        "script": "scripts/prot/tr2aacds4.pl",
        "min_aa": 40,
        "p_heterozygosity": 0,
        "stranded_rna": "yes",
        "aconsensus": True,
        "chunks": 0,
        "chunk_target": 150000,
        "max_merge_levels": 1,
        "chunk_carry_forward": "nonredundant",
        "small_aa_max": None,
        "keep_small_orfs": False,
        "no_utr_orf": False,
        "extra_args": [],
    },
    "tracks": {
        "flnc_prefix": "flnc",
        "dedup_identical": True,
        "cluster": {"enabled": True, "singletons": True, "threads": None,
                    "sort_threads": None, "extra_args": []},
        "keep_raw_flnc_track": False,
        "include_reference_cds": False,
        "ref_cds_prefix": "refcds",
        "drop_reference_cds_when_supported": True,
    },
    "transcriptome": {"include_alternates": True, "min_cds_len": 0,
                      "add_unrecovered_reference": True,
                      "include_cull": False},
    "mapping": {"preset": "splice:hq", "strand_mode": "f", "extra_args": [], "secondary": False},
    "classify": {
        "min_query_cov": 0.50,
        "min_identity": 0.85,
        "min_cds_overlap_bp": 30,
        "min_cds_overlap_frac": 0.0,
        "propagate_ref_by_locus": True,
        "emit_genome_gff": True,
        "quantify_flnc": True,
    },
    "output": {"gff_source": "evigene_pipe", "emit_utr": True, "id_prefix": None},
}


def _deep_merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def _as_list(v) -> list:
    if v is None:
        return []
    return list(v) if isinstance(v, (list, tuple)) else [v]


class Config:
    def __init__(self, data: dict, path: Path | None = None):
        self.path = path
        self.d = data
        self.sample: str = data["sample"]
        # evigene derives file names and public IDs from the input basename; a
        # name containing dots (e.g. "Em2.tr.fa") breaks its stage-2 file lookup.
        self.evg_name: str = re.sub(r"[^A-Za-z0-9]", "", str(self.sample)) or "evg"
        self.outdir: Path = Path(data["outdir"]).expanduser().resolve()
        # The SLURM allocation wins over the config. Without this a job submitted with
        # --cpus-per-task=4 still ran tools with -j/threads from resources.threads,
        # so lowering the sbatch CPU count changed nothing about what the tool did -
        # and per-thread buffers are exactly what makes a memory-hungry step OOM.
        self.threads: int = int(data["resources"]["threads"])
        env_cpus = os.environ.get("SLURM_CPUS_PER_TASK", "")
        if env_cpus.isdigit() and 0 < int(env_cpus) < self.threads:
            LOG.warning("resources.threads=%d exceeds the allocation "
                        "(SLURM_CPUS_PER_TASK=%s) - using %s",
                        self.threads, env_cpus, env_cpus)
            self.threads = int(env_cpus)
        self.max_mem_mb: int = int(data["resources"]["max_mem_mb"])
        env_mem = os.environ.get("SLURM_MEM_PER_NODE", "")
        if not env_mem.isdigit():
            per_cpu = os.environ.get("SLURM_MEM_PER_CPU", "")
            env_mem = str(int(per_cpu) * self.threads) if per_cpu.isdigit() else ""
        if env_mem.isdigit() and 0 < int(env_mem) < self.max_mem_mb:
            LOG.warning("resources.max_mem_mb=%d exceeds the allocation (%s MB) - using %s",
                        self.max_mem_mb, env_mem, env_mem)
            self.max_mem_mb = int(env_mem)
        tmp = data["resources"].get("tmpdir")
        self.tmpdir: Path = Path(tmp).expanduser() if tmp else self.outdir / "tmp"

    # convenience accessors -------------------------------------------------
    @property
    def genome(self) -> Path:
        return Path(self.d["inputs"]["genome"]).expanduser().resolve()

    @property
    def annotation(self) -> Path:
        return Path(self.d["inputs"]["annotation"]).expanduser().resolve()

    @property
    def flnc_bam(self) -> list[Path]:
        return [Path(p).expanduser().resolve() for p in _as_list(self.d["inputs"].get("flnc_bam"))]

    @property
    def flnc_fasta(self) -> list[Path]:
        return [Path(p).expanduser().resolve() for p in _as_list(self.d["inputs"].get("flnc_fasta"))]

    @property
    def evigene_home(self) -> Path:
        return Path(self.d["evigene"]["home"]).expanduser().resolve()

    @property
    def evigene_script(self) -> Path:
        return self.evigene_home / self.d["evigene"]["script"]

    def dir(self, *parts) -> Path:
        p = self.outdir.joinpath(*parts)
        p.mkdir(parents=True, exist_ok=True)
        return p

    def file(self, *parts) -> Path:
        p = self.outdir.joinpath(*parts)
        p.parent.mkdir(parents=True, exist_ok=True)
        return p


def load_config(path: str | Path) -> Config:
    path = Path(path).expanduser().resolve()
    if not path.exists():
        die(f"config not found: {path}")
    with open(path) as fh:
        user = yaml.safe_load(fh) or {}
    data = _deep_merge(DEFAULTS, user)

    errs: list[str] = []
    if not data.get("sample"):
        errs.append("sample: is required")
    if not data.get("outdir"):
        errs.append("outdir: is required")
    for key in ("genome", "annotation"):
        v = data["inputs"].get(key)
        if not v:
            errs.append(f"inputs.{key}: is required")
        elif not Path(v).expanduser().exists():
            errs.append(f"inputs.{key}: no such file: {v}")
    reads = _as_list(data["inputs"].get("flnc_bam")) + _as_list(data["inputs"].get("flnc_fasta"))
    if not reads:
        errs.append("inputs: give at least one of flnc_bam / flnc_fasta")
    for r in reads:
        if not Path(r).expanduser().exists():
            errs.append(f"inputs: no such read file: {r}")
    home = data["evigene"].get("home")
    if not home:
        errs.append("evigene.home: is required (run install_evigene.sh)")
    else:
        script = Path(home).expanduser() / data["evigene"]["script"]
        if not script.exists():
            errs.append(f"evigene.script: not found: {script}")
    ch = data["evigene"].get("chunks")
    if not (str(ch).lower() == "auto" or str(ch).isdigit()):
        errs.append("evigene.chunks must be an integer or \"auto\"")
    if data["tracks"].get("cluster", {}).get("enabled") and not _as_list(
            data["inputs"].get("flnc_bam")):
        errs.append("tracks.cluster.enabled requires inputs.flnc_bam "
                    "(isoseq cluster2 cannot read fasta); set it to false for fasta input")
    if data["mapping"]["strand_mode"] not in ("f", "b", "n"):
        errs.append("mapping.strand_mode must be one of f, b, n")
    if str(data["evigene"]["stranded_rna"]) not in ("auto", "yes", "no"):
        errs.append("evigene.stranded_rna must be auto, yes or no")

    if errs:
        die("config errors in %s:\n  - %s" % (path, "\n  - ".join(errs)))

    LOG.debug("config loaded from %s", path)
    return Config(data, path)
