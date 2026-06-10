#!/usr/bin/env bash
# run.sh — End-to-end EIC FFS study campaign (ANALYSIS_DESIGN.md)
#
#   Stage 1  generate_events.py   Pythia8 NC-DIS, 3 beam configs +
#                                 hadronization variations + MPI check
#   Stage 2  analyze_jets.py      gamma*p-frame jets -> per-jet tables
#                                 (particle level + smeared reco level)
#   Stage 3  make_results.py      profiles, H0 test, collapse metric,
#                                 luminosity projections -> results.json
#   Stage 4  make_paper_plots.py  the four PRL figures
#
# Usage:
#   ./run.sh                 # full campaign (~3M events, ~30 min on 4 cores)
#   SCALE=0.1 ./run.sh       # 10% statistics quick run
#   ./run.sh --skip-gen      # re-analyze existing events
#   ./run.sh --skip-ana      # re-run stats/plots only

set -euo pipefail

SCALE="${SCALE:-1.0}"
DATA_DIR="${DATA_DIR:-data}"
PLOT_DIR="${PLOT_DIR:-plots}"
NPROC="${NPROC:-$(nproc)}"

SKIP_GEN=false
SKIP_ANA=false
for arg in "$@"; do
    case $arg in
        --skip-gen) SKIP_GEN=true ;;
        --skip-ana) SKIP_ANA=true ;;
    esac
done

nev() { python3 -c "print(max(1000, int($1 * $SCALE)))"; }

mkdir -p "${DATA_DIR}" "${PLOT_DIR}" logs results

# ── Stage 1: generation campaign ────────────────────────────────────────────
if [ "$SKIP_GEN" = false ]; then
    echo "── Stage 1: Pythia8 generation campaign ──"
    JOBS=$(cat <<EOF
--config 5x41 --variation baseline --n-events $(nev 500000) --seed 101
--config 10x100 --variation baseline --n-events $(nev 500000) --seed 102
--config 18x275 --variation baseline --n-events $(nev 500000) --seed 103
--config 5x41 --variation lund-soft --n-events $(nev 200000) --seed 104
--config 5x41 --variation lund-hard --n-events $(nev 200000) --seed 105
--config 10x100 --variation lund-soft --n-events $(nev 200000) --seed 106
--config 10x100 --variation lund-hard --n-events $(nev 200000) --seed 107
--config 18x275 --variation lund-soft --n-events $(nev 200000) --seed 108
--config 18x275 --variation lund-hard --n-events $(nev 200000) --seed 109
--config 10x100 --variation baseline --mpi --n-events $(nev 150000) --seed 110
EOF
)
    echo "$JOBS" | xargs -P "$NPROC" -I{} sh -c \
        'tag=$(echo "{}" | tr -cd "a-zA-Z0-9" | tail -c 40); \
         python3 generate_events.py {} --quiet > logs/gen_$tag.log 2>&1'
fi

# ── Stage 1b: Herwig leg (cluster hadronization), if installed ──────────────
if [ "$SKIP_GEN" = false ]; then
    if command -v Herwig > /dev/null 2>&1 || [ -x "${HEP_PREFIX:-/opt/hep}/bin/Herwig" ]; then
        echo "── Stage 1b: Herwig 7 campaign ──"
        ./run_herwig.sh "$(nev 300000)"
    else
        echo "── Stage 1b: Herwig not found; skipping cluster-hadronization leg ──"
    fi
fi

# ── Stage 2: jet analysis (truth + smeared reco) ────────────────────────────
if [ "$SKIP_ANA" = false ]; then
    echo "── Stage 2: jet analysis ──"
    { ls ${DATA_DIR}/events_*.parquet;
      for c in 5x41 10x100 18x275; do
          echo "${DATA_DIR}/events_${c}_baseline.parquet --smear"; done; } |
    xargs -P "$NPROC" -I{} sh -c \
        'tag=$(echo "{}" | tr -cd "a-zA-Z0-9" | tail -c 44); \
         python3 analyze_jets.py {} > logs/ana_$tag.log 2>&1'
fi

# ── Stage 3 + 4: statistics and figures ─────────────────────────────────────
echo "── Stage 3: statistics ──"
python3 make_results.py --datadir "${DATA_DIR}" --output results/results.json

echo "── Stage 4: figures ──"
python3 make_paper_plots.py --results results/results.json --outdir "${PLOT_DIR}"

echo ""
echo "Campaign complete:  results/results.json  +  ${PLOT_DIR}/fig*.pdf"
