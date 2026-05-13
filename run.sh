#!/usr/bin/env bash
# run.sh — End-to-end EIC FFS study pipeline
#
# Usage:
#   ./run.sh                        # full run (200k events)
#   N_EVENTS=50000 ./run.sh         # quick test run
#   ./run.sh --skip-gen             # skip generation, re-analyze existing events
#   ./run.sh --skip-gen --skip-ana  # re-plot only

set -euo pipefail

# ── Configuration ──────────────────────────────────────────────────────────
N_EVENTS="${N_EVENTS:-200000}"
ELECTRON_ENERGY="${ELECTRON_ENERGY:-10.0}"
PROTON_ENERGY="${PROTON_ENERGY:-100.0}"
Q2MIN="${Q2MIN:-1.0}"
DATA_DIR="${DATA_DIR:-data}"
PLOT_DIR="${PLOT_DIR:-plots}"
EVENTS_FILE="${DATA_DIR}/events.parquet"
HISTO_FILE="${DATA_DIR}/histograms.root"

SKIP_GEN=false
SKIP_ANA=false
USE_FASTJET=""

for arg in "$@"; do
    case $arg in
        --skip-gen)  SKIP_GEN=true ;;
        --skip-ana)  SKIP_ANA=true ;;
        --fastjet)   USE_FASTJET="--use-fastjet" ;;
    esac
done

echo "================================================================"
echo "  EIC FFS Study Pipeline"
echo "  Beam: ${ELECTRON_ENERGY} GeV e  ×  ${PROTON_ENERGY} GeV p"
echo "  Events: ${N_EVENTS}"
echo "================================================================"
echo ""

mkdir -p "${DATA_DIR}" "${PLOT_DIR}"

# ── Step 1: Event generation ───────────────────────────────────────────────
if [ "$SKIP_GEN" = false ]; then
    echo "── Step 1: Generating ${N_EVENTS} DIS events with Pythia8 ─────────"
    python generate_events.py \
        --n-events       "${N_EVENTS}" \
        --electron-energy "${ELECTRON_ENERGY}" \
        --proton-energy   "${PROTON_ENERGY}" \
        --Q2min           "${Q2MIN}" \
        --output          "${EVENTS_FILE}" \
        --quiet
    echo ""
else
    echo "── Step 1: Skipping event generation (using ${EVENTS_FILE}) ────────"
fi

if [ ! -f "${EVENTS_FILE}" ]; then
    echo "ERROR: ${EVENTS_FILE} not found. Run without --skip-gen first." >&2
    exit 1
fi

# ── Step 2: Analysis ───────────────────────────────────────────────────────
if [ "$SKIP_ANA" = false ]; then
    echo "── Step 2: Analyzing events and filling histograms ─────────────────"
    python analyze_events.py \
        "${EVENTS_FILE}" \
        --output "${HISTO_FILE}" \
        ${USE_FASTJET:+"${USE_FASTJET}"}
    echo ""
else
    echo "── Step 2: Skipping analysis (using ${HISTO_FILE}) ─────────────────"
fi

if [ ! -f "${HISTO_FILE}" ]; then
    echo "ERROR: ${HISTO_FILE} not found." >&2
    exit 1
fi

# ── Step 3: Plotting ───────────────────────────────────────────────────────
echo "── Step 3: Making publication-quality plots ─────────────────────────"
python make_plots.py \
    "${HISTO_FILE}" \
    --outdir "${PLOT_DIR}"

echo ""
echo "================================================================"
echo "  Pipeline complete!"
echo "  Histograms: ${HISTO_FILE}"
echo "  Plots:      ${PLOT_DIR}/"
echo "================================================================"
