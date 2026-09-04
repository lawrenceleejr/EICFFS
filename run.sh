#!/usr/bin/env bash
# run.sh — End-to-end EIC FFS study pipeline
#
# Usage:
#   ./run.sh                          # full run: N_SEEDS × N_EVENTS events, analysis, figures
#   N_EVENTS=50000 N_SEEDS=1 ./run.sh # quick test run
#   ./run.sh --skip-gen               # re-analyse existing events
#   ./run.sh --skip-gen --skip-ana    # re-plot only
#   ./run.sh --legacy-plots           # also produce the original make_plots.py figures
#
# Generation runs N_SEEDS Pythia jobs in parallel (N_PARALLEL at a time), each
# saving N_EVENTS events with W > WMIN, to data/events_<seed>.parquet.

set -euo pipefail

# ── Configuration ──────────────────────────────────────────────────────────
N_EVENTS="${N_EVENTS:-300000}"
N_SEEDS="${N_SEEDS:-8}"
N_PARALLEL="${N_PARALLEL:-4}"
ELECTRON_ENERGY="${ELECTRON_ENERGY:-10.0}"
PROTON_ENERGY="${PROTON_ENERGY:-100.0}"
Q2MIN="${Q2MIN:-1.0}"
Q2MAX="${Q2MAX:-1000.0}"
WMIN="${WMIN:-10.0}"
DATA_DIR="${DATA_DIR:-data}"
FIG_DIR="${FIG_DIR:-figures}"
PLOT_DIR="${PLOT_DIR:-plots}"
ANALYSIS_FILE="${DATA_DIR}/analysis.root"

SKIP_GEN=false
SKIP_ANA=false
LEGACY=false
FASTJET_FLAG=""

for arg in "$@"; do
    case $arg in
        --skip-gen)     SKIP_GEN=true ;;
        --skip-ana)     SKIP_ANA=true ;;
        --legacy-plots) LEGACY=true ;;
        --no-fastjet)   FASTJET_FLAG="--no-fastjet" ;;
    esac
done

echo "================================================================"
echo "  EIC FFS Study Pipeline"
echo "  Beam: ${ELECTRON_ENERGY} GeV e  ×  ${PROTON_ENERGY} GeV p"
echo "  Events: ${N_SEEDS} × ${N_EVENTS}  (W > ${WMIN} GeV)"
echo "================================================================"
echo ""

mkdir -p "${DATA_DIR}" "${FIG_DIR}"

# ── Step 1: Event generation ───────────────────────────────────────────────
if [ "$SKIP_GEN" = false ]; then
    echo "── Step 1: Generating events with Pythia8 (${N_PARALLEL} parallel jobs) ──"
    seq 1 "${N_SEEDS}" | xargs -P "${N_PARALLEL}" -I{} sh -c \
        "python3 generate_events.py \
            --n-events ${N_EVENTS} --seed {} --Wmin ${WMIN} \
            --electron-energy ${ELECTRON_ENERGY} --proton-energy ${PROTON_ENERGY} \
            --Q2min ${Q2MIN} --Q2max ${Q2MAX} \
            --output ${DATA_DIR}/events_{}.parquet --quiet \
            > ${DATA_DIR}/gen_{}.log 2>&1 && echo '  seed {} done'"
    echo ""
else
    echo "── Step 1: Skipping event generation ────────────────────────────────"
fi

if ! ls "${DATA_DIR}"/events_*.parquet >/dev/null 2>&1; then
    echo "ERROR: no ${DATA_DIR}/events_*.parquet found. Run without --skip-gen first." >&2
    exit 1
fi

# ── Step 2: Analysis ───────────────────────────────────────────────────────
if [ "$SKIP_ANA" = false ]; then
    echo "── Step 2: Jet finding, frame boosts, n90 → ${ANALYSIS_FILE} ─────────"
    python3 analyze_events.py "${DATA_DIR}"/events_*.parquet \
        --output "${ANALYSIS_FILE}" ${FASTJET_FLAG}
    echo ""
else
    echo "── Step 2: Skipping analysis (using ${ANALYSIS_FILE}) ───────────────"
fi

# ── Step 3: Figures ────────────────────────────────────────────────────────
echo "── Step 3: Figures → ${FIG_DIR}/ ────────────────────────────────────"
FIRST_EVENTS=$(ls "${DATA_DIR}"/events_*.parquet | head -1)
python3 make_figures.py "${ANALYSIS_FILE}" --events "${FIRST_EVENTS}" --outdir "${FIG_DIR}"

if [ "$LEGACY" = true ]; then
    mkdir -p "${PLOT_DIR}"
    python3 make_plots.py "${ANALYSIS_FILE}" --outdir "${PLOT_DIR}"
fi

echo ""
echo "================================================================"
echo "  Pipeline complete!"
echo "  Analysis: ${ANALYSIS_FILE}"
echo "  Figures:  ${FIG_DIR}/"
echo "================================================================"
