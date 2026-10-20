#!/usr/bin/env bash
###############################################################################
##  run_local.sh -- runs the whole simulation study (studies 1 and 2)
##  locally in parallel, inside the pinned Nix environment (see gen-env.R).
##
##  Usage (from the repository root, inside `nix-shell --pure` or via direnv):
##    JOBS=12 bash simulation_study/run_local.sh          # run all cells
##    bash simulation_study/run_local.sh --merge-only     # only merge results
##
##  Same cells as the run_study_{1,2}.sh SLURM array jobs (both use cells.sh and
##  run_cell.sh). Each cell (one call of run_simulation.R = 100 replicates)
##  writes its own CSV under
##  simulation_study/results/cells/, so parallel jobs never append to the
##  same file; finished cells are skipped when the script is restarted.
##  At the end the cell files are merged into the per-study result files that
##  the plotting scripts read.
###############################################################################
set -euo pipefail
cd "$(dirname "$0")/.."

JOBS="${JOBS:-12}"
CELLS="simulation_study/results/cells"
# one core per cell; avoid BLAS oversubscription
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [[ "${1:-}" != "--merge-only" ]]; then
  mkdir -p "$CELLS"
  # cells.sh lists the design (shared with the SLURM array jobs); run_cell.sh
  # runs one cell and skips it if it is already finished
  bash simulation_study/cells.sh > "$CELLS/cells.txt"
  echo "$(wc -l < "$CELLS/cells.txt") cells, $JOBS parallel jobs ($(date))"
  xargs -P "$JOBS" -L 1 bash simulation_study/run_cell.sh < "$CELLS/cells.txt" |
    grep -v "^already done" || true
  echo "all cells finished ($(date))"
fi

# check and merge the cell files into the per-study result files
bash simulation_study/merge_results.sh
