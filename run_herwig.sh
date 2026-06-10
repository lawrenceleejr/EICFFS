#!/usr/bin/env bash
# run_herwig.sh — Herwig 7 leg of the FFS campaign (cluster hadronization).
#
# For each EIC beam configuration: fill the run-card template, integrate,
# generate N events with HepMC output, parse the generated cross section,
# convert to the EICFFS Parquet schema, and remove the bulky HepMC file.
#
# Usage:
#   ./run_herwig.sh [N_EVENTS] [config ...]      # default 300000, all configs
#
# Requires Herwig in PATH (default /opt/hep) and pyhepmc (pip).

set -euo pipefail

HEP_PREFIX="${HEP_PREFIX:-/opt/hep}"
export PATH="$HEP_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$HEP_PREFIX/lib:$HEP_PREFIX/lib64:${LD_LIBRARY_PATH:-}"
export LHAPDF_DATA_PATH="$HEP_PREFIX/share/LHAPDF:${LHAPDF_DATA_PATH:-}"

N_EVENTS="${1:-300000}"
shift || true
CONFIGS=("$@")
[ ${#CONFIGS[@]} -eq 0 ] && CONFIGS=(5x41 10x100 18x275)

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKBASE="${WORKBASE:-/tmp/herwig_runs}"
mkdir -p "$WORKBASE" "$REPO_DIR/data" "$REPO_DIR/logs"

beam_energies() {
    case "$1" in
        5x41)   echo "5 41" ;;
        10x100) echo "10 100" ;;
        18x275) echo "18 275" ;;
        *) echo "unknown config $1" >&2; exit 1 ;;
    esac
}

seed_for() {
    case "$1" in
        5x41) echo 911 ;; 10x100) echo 912 ;; 18x275) echo 913 ;;
    esac
}

run_one() {
    local config="$1"
    read -r EE EP <<< "$(beam_energies "$config")"
    local wd="$WORKBASE/$config"
    local run="EIC_${config}"
    local hepmc="$wd/$run.hepmc"
    mkdir -p "$wd"

    sed -e "s/@EE@/$EE/" -e "s/@EP@/$EP/" \
        -e "s/@Q2MIN@/25/" -e "s/@Q2MAX@/1000/" \
        -e "s#@HEPMC@#$hepmc#" -e "s/@RUN@/$run/" \
        "$REPO_DIR/herwig/EIC_DIS.in.template" > "$wd/$run.in"

    echo "[$config] Herwig read"
    (cd "$wd" && Herwig read "$run.in" > read.log 2>&1)

    echo "[$config] Herwig run: $N_EVENTS events"
    local seed
    seed="$(seed_for "$config")"
    (cd "$wd" && Herwig run "$run.run" -N "$N_EVENTS" \
        -s "$seed" -d 0 > run.log 2>&1)

    # Parse generated cross section from the ThePEG .out file -> fb
    local sigma_fb
    sigma_fb=$(python3 - "$wd/$run-S$seed.out" <<'PYEOF'
import re, sys
txt = open(sys.argv[1]).read()
# Herwig/ThePEG format: "Total (from attempted events): ... 12.4(1)e+00"
# in the units given in the table header (nb).
m = re.search(r"Total \(from attempted events\).*?([0-9.]+)\(\d+\)e([+-]\d+)",
              txt)
if not m:
    sys.exit("could not parse cross section from " + sys.argv[1])
val = float(m.group(1)) * 10.0**int(m.group(2))     # nb
print(val * 1e6)                                     # fb
PYEOF
)
    echo "[$config] sigma = $sigma_fb fb"

    python3 "$REPO_DIR/convert_hepmc.py" "$hepmc" --config "$config" \
        --output "$REPO_DIR/data/events_${config}_herwig.parquet" \
        --sigma-fb "$sigma_fb" > "$REPO_DIR/logs/convert_${config}.log" 2>&1
    rm -f "$hepmc"
    echo "[$config] DONE"
}

for config in "${CONFIGS[@]}"; do
    run_one "$config" &
done
wait
echo "HERWIG_CAMPAIGN_DONE"
