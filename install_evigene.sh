#!/usr/bin/env bash
# Install EvidentialGene (not available on conda) and check tr2aacds dependencies.
#
#   bash install_evigene.sh [INSTALL_DIR]
#
# Default INSTALL_DIR is ./software . Afterwards set evigene_home in your config
# to the printed path, or `export evigene=<path>`.
set -euo pipefail

INSTALL_DIR="${1:-$(pwd)/software}"
TARBALL_URL="${EVIGENE_URL:-https://downloads.sourceforge.net/project/evidentialgene/evigene23jul15.tar}"
TARBALL_NAME="$(basename "${TARBALL_URL%%\?*}")"

mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

if [[ -d evigene/scripts/prot ]]; then
  echo "[install_evigene] evigene already present at $INSTALL_DIR/evigene"
else
  echo "[install_evigene] downloading $TARBALL_URL"
  curl -fsSL -o "$TARBALL_NAME" "$TARBALL_URL"
  echo "[install_evigene] unpacking $TARBALL_NAME"
  tar -xf "$TARBALL_NAME"
  rm -f "$TARBALL_NAME"
fi

EVIGENE_HOME="$INSTALL_DIR/evigene"
[[ -d "$EVIGENE_HOME" ]] || { echo "ERROR: $EVIGENE_HOME missing after unpack" >&2; exit 1; }
chmod -R u+rx "$EVIGENE_HOME/scripts" 2>/dev/null || true

# --- make evigene use the perl on PATH -------------------------------------
# tr2aacds runs its helper scripts as executables (e.g. `$EVIGENES/cdna_bestorf.pl
# -cdna ...`), never as `perl <script>`, so the shebang decides the interpreter.
# 265 of them hardcode #!/usr/bin/perl, which bypasses the conda environment and
# lands on the system perl. On RHEL/Rocky the system perl is modular and core
# modules like FindBin live in a separate RPM that is often absent, giving
# "Can't locate FindBin.pm in @INC" partway into a long job. Rewriting to
# /usr/bin/env perl makes PATH decide, so the conda perl is used throughout.
# match a direct interpreter path ending in /perl; '#!/usr/bin/env perl' has a
# space before perl and so is correctly excluded
# grep exits 1 when nothing matches, which under `set -e -o pipefail` would abort
# the script on every re-run once the shebangs are already fixed
n_rw=$({ grep -rlE '^#!(/[^ ]*)*/perl[[:space:]]*(-w)?[[:space:]]*$' "$EVIGENE_HOME/scripts" \
         --include='*.pl' 2>/dev/null || true; } | wc -l)
if [[ "$n_rw" -gt 0 ]]; then
  find "$EVIGENE_HOME/scripts" -name '*.pl' -exec \
    sed -i '1s|^#!.*[/ ]perl[[:space:]]*\(-w\)\?[[:space:]]*$|#!/usr/bin/env perl|' {} +
  echo "[install_evigene] rewrote $n_rw perl shebangs to /usr/bin/env perl"
else
  echo "[install_evigene] perl shebangs already use /usr/bin/env"
fi

echo
echo "[install_evigene] tr2aacds entry points found:"
ls "$EVIGENE_HOME"/scripts/prot/tr2aacds*.pl | sed 's/^/    /'

echo
echo "[install_evigene] dependency check:"
missing=0
for tool in perl fastanrdb cd-hit cd-hit-est blastn makeblastdb minimap2 samtools gffread; do
  if command -v "$tool" >/dev/null 2>&1; then
    printf '    %-14s OK  %s\n' "$tool" "$(command -v "$tool")"
  else
    printf '    %-14s MISSING\n' "$tool"
    missing=1
  fi
done
python - <<'PY' || missing=1
import importlib.util as iu
import sys
mods = ("pysam", "Bio", "pandas", "yaml")
bad = [m for m in mods if iu.find_spec(m) is None]
for m in mods:
    print(f"    {m:<14} {'OK' if m not in bad else 'MISSING'}")
sys.exit(1 if bad else 0)
PY

echo
echo "[install_evigene] running an evigene helper script for real:"
bestorf="$EVIGENE_HOME/scripts/cdna_bestorf.pl"
if out=$("$bestorf" 2>&1); then :; fi
if grep -qE "Can't locate|BEGIN failed|command not found" <<<"$out"; then
  echo "    FAILED - the helper scripts cannot start:"
  sed 's/^/      /' <<<"$out" | head -5
  echo "    The perl running them is: $(head -1 "$bestorf")"
  missing=1
elif grep -q "usage:" <<<"$out"; then
  echo "    OK  ($("$bestorf" -showperl 2>/dev/null || command -v perl))"
else
  echo "    unexpected output from $bestorf:"
  sed 's/^/      /' <<<"$out" | head -5
fi

echo
if [[ $missing -ne 0 ]]; then
  echo "[install_evigene] Some dependencies are missing. Activate the conda env first:"
  echo "    conda env create -f environment.yml && conda activate evigene_pipe"
  exit 1
fi

echo "[install_evigene] done."
echo "    evigene_home: $EVIGENE_HOME"
echo "    export evigene=$EVIGENE_HOME"
