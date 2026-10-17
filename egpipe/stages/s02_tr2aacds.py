"""Stage 02 - run EvidentialGene tr2aacds on the merged input fasta.

tr2aacds writes its output folders (okayset/, dropset/, inputset/, tmpfiles/)
into the current working directory, so it is run inside 02_evigene/.
"""
from __future__ import annotations

import re
import shutil
from collections import Counter
from pathlib import Path

from ..utils import LOG, die, iter_fasta, require_tools, run, write_fasta

STAGE = "02_tr2aacds"


def build_command(cfg, input_fa: Path, species: str | None = None) -> list[str]:
    ev = cfg.d["evigene"]
    cmd: list[str] = [
        "perl", str(cfg.evigene_script),
        "-cdnaseq", input_fa.name,
        "-NCPU", str(cfg.threads),
        "-MAXMEM", str(cfg.max_mem_mb),
        "-MINAA", str(ev["min_aa"]),
        "-logfile",
        "-tidyup",
    ]
    if int(ev.get("p_heterozygosity", 0)) > 0:
        cmd += ["-pHeterozygosity", str(int(ev["p_heterozygosity"]))]
    if str(ev.get("stranded_rna", "auto")) != "auto":
        cmd += ["-strandedrna", str(ev["stranded_rna"])]
    cmd += ["-aconsensus"] if ev.get("aconsensus", True) else ["-aconsensus=no"]
    if ev.get("small_aa_max") is not None:
        cmd += ["-SMALLAAMAX", str(int(ev["small_aa_max"]))]
    if ev.get("keep_small_orfs"):
        cmd += ["-smallclass"]
    if ev.get("no_utr_orf"):
        cmd += ["-noutrorf"]
    cmd += ["-species", str(species or cfg.evg_name)]
    cmd += [str(a) for a in (ev.get("extra_args") or [])]
    return cmd


def check_evigene_perl(cfg) -> None:
    """Actually run an evigene helper script before committing to a long job.

    tr2aacds executes its helpers as executables, so their `#!` line - not PATH,
    and not the perl this pipeline was started with - picks the interpreter. The
    stock package hardcodes #!/usr/bin/perl, which on RHEL/Rocky finds a modular
    system perl whose core modules (FindBin, File::Basename) are separate RPMs
    that are frequently not installed. That surfaces as
    "Can't locate FindBin.pm in @INC" from a helper, often well into the run.
    """
    import subprocess

    helper = cfg.evigene_home / "scripts" / "cdna_bestorf.pl"
    if not helper.exists():
        die(f"evigene helper not found: {helper}")
    shebang = helper.read_text(errors="replace").split("\n", 1)[0].strip()
    try:
        proc = subprocess.run([str(helper)], capture_output=True, text=True, timeout=60)
    except PermissionError:
        die(f"{helper} is not executable - run: chmod -R u+rx {cfg.evigene_home}/scripts")
    except OSError as e:
        die(f"cannot execute {helper}: {e}\n  its shebang is: {shebang}")
    out = (proc.stdout or "") + (proc.stderr or "")

    if "Can't locate" in out or "BEGIN failed" in out:
        first = next((l for l in out.splitlines() if "Can't locate" in l), out.splitlines()[0])
        die(
            f"evigene helper scripts cannot start.\n"
            f"  {first.strip()}\n"
            f"  {helper.name} runs under: {shebang}\n"
            f"  tr2aacds calls its helpers as executables, so that shebang - not the\n"
            f"  conda environment - chooses the perl. Point them at the perl on PATH:\n"
            f"    find {cfg.evigene_home}/scripts -name '*.pl' -exec \\\n"
            f"      sed -i '1s|^#!.*[/ ]perl.*$|#!/usr/bin/env perl|' {{}} +\n"
            f"  (install_evigene.sh now does this automatically.)"
        )
    if "usage:" not in out:
        LOG.warning("unexpected output from %s (continuing):\n    %s",
                    helper.name, "\n    ".join(out.splitlines()[:4]))
    else:
        LOG.info("evigene helpers OK (%s runs under %s)", helper.name, shebang)


def _evigene_env(cfg) -> dict:
    cfg.tmpdir.mkdir(parents=True, exist_ok=True)
    return {
        "EVIGENES": str(cfg.evigene_home / "scripts"),
        "evigene": str(cfg.evigene_home),
        "TMPDIR": str(cfg.tmpdir),
    }


def _log_offsets(work: Path) -> dict:
    """Byte size of each evigene log before a run.

    evigene opens its log with '>>', so a directory that has been run before
    carries every earlier attempt's messages - including the FATAL from a run
    that died at step 0. Scanning the whole file reports those stale errors
    against a run that actually succeeded, so only the bytes this run appended
    are examined.
    """
    return {p: p.stat().st_size for p in work.glob("*.log") if p.is_file()}


def _check_outputs(work: Path, tag: str, log_offsets: dict | None = None) -> Path:
    okayset = work / "okayset"
    fatal = []
    for lg in sorted(work.glob("*.log")):
        start = (log_offsets or {}).get(lg, 0)
        try:
            with open(lg, "rb") as fh:
                if start and start <= lg.stat().st_size:
                    fh.seek(start)
                text = fh.read().decode("utf-8", errors="replace")
        except OSError:
            continue
        for line in text.splitlines():
            if "FATAL" in line or line.startswith("#egr: ERR"):
                fatal.append(f"{lg.name}: {line.strip()}")
    usable = okayset.is_dir() and (any(okayset.glob("*.mrna")) or any(okayset.glob("*.okay.*")))
    if not usable:
        msg = f"tr2aacds produced no usable okayset in {work} ({tag})"
        if fatal:
            msg += "\n  evigene reported:\n    " + "\n    ".join(fatal[:10])
        die(msg)
    if fatal:
        LOG.warning("evigene reported problems in %s (outputs still present):\n    %s",
                    tag, "\n    ".join(fatal[:10]))
    return work


def _split_fasta(input_fa: Path, outdir: Path, n_chunks: int, basename: str) -> list[Path]:
    """Round-robin split so chunks are balanced in count and length."""
    handles, paths = [], []
    for k in range(n_chunks):
        d = outdir / f"chunk{k + 1:02d}"
        d.mkdir(parents=True, exist_ok=True)
        p = d / f"{basename}c{k + 1}.fa"
        paths.append(p)
        handles.append(open(p, "w"))
    try:
        for i, (name, _desc, seq) in enumerate(iter_fasta(input_fa)):
            write_fasta(handles[i % n_chunks], name, seq)
    finally:
        for h in handles:
            h.close()
    for p in paths:
        LOG.info("  chunk %s: %d sequences", p.parent.name,
                 sum(1 for line in open(p) if line.startswith(">")))
    return paths


def _pubid_to_oids(evgdir: Path) -> dict[str, str]:
    """chunk public ID -> the original input IDs behind it (from the pubids table)."""
    out: dict[str, str] = {}
    cands = [c for c in list((evgdir / "okayset").glob("*.pubids"))
             + list((evgdir / "publicset").glob("*.pubids")) if not c.name.endswith(".old")]
    if not cands:
        return out
    with open(cands[0]) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            oids = f[8] if len(f) > 8 and f[8] else f[1]
            out[f[0]] = oids
    return out


def _count_fasta(path: Path) -> int:
    return sum(1 for line in open(path) if line.startswith(">"))


def _expand(ids: str, chain: dict[str, str]) -> str:
    """Resolve a comma list of level-N ids down to the original input ids."""
    if not chain or not ids:
        return ids
    out, seen = [], set()
    for i in str(ids).split(","):
        i = i.strip()
        if not i:
            continue
        for x in chain.get(i, i).split(","):
            if x and x not in seen:
                seen.add(x)
                out.append(x)
    return ",".join(out)


# Drop reasons that are safe to discard inside a chunk: the sequence is byte-identical
# to, or fully contained in, another sequence in the same chunk, so nothing is lost by
# not carrying it forward. Every other drop reason (parthi, altmap, smallaa, althi, ...)
# is a *classification* call made without global context - those are carried forward so
# the final combine decides them against the whole set.
DROP_REDUNDANT = ("perfectdup", "perfectfrag")


def _dropset_carry(evgdir: Path):
    """Yield (name, seq) for dropped transcripts that deserve another look."""
    d = evgdir / "dropset"
    if not d.is_dir():
        return
    files = sorted(set(list(d.glob("*.drop.tr")) + list(d.glob("*.drop.tr.fa"))
                       + list(d.glob("*.drop.mrna"))))
    for f in files:
        for name, desc, seq in iter_fasta(f):
            m = re.search(r"evgclass=([A-Za-z0-9]+)", desc or "")
            cls = m.group(1) if m else ""
            if cls.startswith(DROP_REDUNDANT):
                continue
            yield name, seq, cls


def _run_chunked(cfg, input_fa: Path, n_chunks: int) -> Path:
    """Reduce in levels, then combine once the set fits in a single run.

    evigene documents running tr2aacds on subsets and re-running on the okay
    outputs, and its classifier is deterministic with respect to how inputs are
    grouped. A single merge round is not enough at this scale: K chunks each
    keeping 30-70% leaves a combine step several times larger than any one chunk,
    which just moves the bottleneck. So the merge repeats until the surviving set
    is under chunk_target, and only then is combined.
    """
    from .s03_harvest import _find_seqsets

    work = cfg.dir("02_evigene")
    env = _evigene_env(cfg)
    target = int(cfg.d["evigene"].get("chunk_target") or 150000)
    max_levels = int(cfg.d["evigene"].get("max_merge_levels") or 4)

    carry_drops = str(cfg.d["evigene"].get("chunk_carry_forward", "nonredundant")) != "okay"
    drop_classes: Counter = Counter()
    chain: dict[str, str] = {}          # current-level id -> ORIGINAL input ids
    current = input_fa
    n_cur = _count_fasta(current)
    level, k = 1, n_chunks

    while k > 1:
        if level > max_levels:
            LOG.warning("reached max_merge_levels=%d with %d sequences still to combine; "
                        "the final run will be large", max_levels, n_cur)
            break
        leveldir = work / f"level{level}"
        LOG.info("=== level %d: %d sequences into %d chunks ===", level, n_cur, k)
        chunk_fas = _split_fasta(current, leveldir, k, f"{cfg.evg_name}L{level}")

        merged = leveldir / "merged.fa"
        new_chain: dict[str, str] = {}
        n_out = 0
        with open(merged, "w") as out:
            for ci, fa in enumerate(chunk_fas, 1):
                d = fa.parent
                tag = f"level{level}/chunk{ci:02d}"
                offsets = _log_offsets(d)
                if (d / "okayset").is_dir() and not getattr(cfg, "force", False):
                    LOG.info("%s already has an okayset - skipping", tag)
                else:
                    LOG.info("--- tr2aacds %s (%d/%d) ---", tag, ci, k)
                    run(build_command(cfg, fa, species=f"{cfg.evg_name}L{level}c{ci}"),
                        cwd=d, env=env)
                _check_outputs(d, tag, offsets)
                sets = _find_seqsets(d, include_alt=True, include_cull=True)
                if not sets.get("mrna"):
                    die(f"{tag}: no okayset mRNA found")
                pubmap = _pubid_to_oids(d)
                j = 0
                carried = 0
                # Sequences are renamed here rather than reusing evigene's public IDs:
                # -species is truncated to 6 characters when the ID prefix is built, so
                # per-chunk species names still collide and the merge would abort on
                # duplicate IDs.
                for part in sets["mrna"]:
                    for name, _desc, seq in iter_fasta(part):
                        j += 1
                        n_out += 1
                        new = f"L{level}c{ci}t{j:07d}"
                        new_chain[new] = _expand(pubmap.get(name, name), chain)
                        write_fasta(out, new, seq)
                if carry_drops:
                    for name, seq, cls in _dropset_carry(d):
                        j += 1
                        n_out += 1
                        carried += 1
                        drop_classes[cls] += 1
                        new = f"L{level}c{ci}t{j:07d}"
                        new_chain[new] = _expand(pubmap.get(name, name), chain)
                        write_fasta(out, new, seq)
                LOG.info("  %s kept %d transcripts (%d carried from dropset)",
                         tag, j, carried)

        pct = 100.0 * n_out / n_cur if n_cur else 0.0
        LOG.info("=== level %d done: %d -> %d sequences (%.1f%% kept) ===",
                 level, n_cur, n_out, pct)
        if drop_classes:
            LOG.info("    carried forward from dropsets: %s",
                     ", ".join(f"{k}={v}" for k, v in drop_classes.most_common(8)))
            drop_classes.clear()
        chain = new_chain
        if n_out >= n_cur:
            LOG.warning("level %d did not reduce the set - stopping the cascade", level)
            current, n_cur = merged, n_out
            break
        current, n_cur = merged, n_out
        level += 1
        k = max(1, -(-n_cur // target))

    chain_tsv = work / "oid_chain.tsv"
    with open(chain_tsv, "w") as fh:
        fh.write("merge_id\toriginal_oids\n")
        for kk, v in sorted(chain.items()):
            fh.write(f"{kk}\t{v}\n")
    LOG.info("ID chain: %d entries -> %s", len(chain), chain_tsv.name)

    final_dir = work / "final"
    final_dir.mkdir(parents=True, exist_ok=True)
    final_fa = final_dir / f"{cfg.evg_name}.fa"
    if current != final_fa:
        shutil.copyfile(current, final_fa)
    LOG.info("=== final combine: %d sequences ===", n_cur)
    offsets = _log_offsets(final_dir)
    if not (final_dir / "okayset").is_dir() or getattr(cfg, "force", False):
        run(build_command(cfg, final_fa), cwd=final_dir, env=env)
    return _check_outputs(final_dir, "final", offsets)


def _resolve_chunks(cfg, input_fa: Path) -> int:
    spec = cfg.d["evigene"].get("chunks") or 0
    if str(spec).lower() == "auto":
        target = int(cfg.d["evigene"].get("chunk_target") or 150000)
        n = sum(1 for line in open(input_fa) if line.startswith(">"))
        k = max(1, -(-n // target))
        LOG.info("chunks=auto: %d sequences / target %d -> %d chunk(s)", n, target, k)
        return k
    return int(spec)


def run_stage(cfg, input_fa: Path) -> Path:
    require_tools(["perl", "fastanrdb", "cd-hit", "cd-hit-est", "blastn", "makeblastdb"])
    check_evigene_perl(cfg)
    n_chunks = _resolve_chunks(cfg, input_fa)
    if n_chunks > 1:
        return _run_chunked(cfg, input_fa, n_chunks)

    work = cfg.dir("02_evigene")
    local = work / f"{cfg.evg_name}.fa"
    if not local.exists() or local.stat().st_size != input_fa.stat().st_size:
        shutil.copyfile(input_fa, local)

    LOG.info("tr2aacds working dir: %s", work)
    offsets = _log_offsets(work)
    run(build_command(cfg, local), cwd=work, env=_evigene_env(cfg))

    return _check_outputs(work, "single-pass", offsets)
