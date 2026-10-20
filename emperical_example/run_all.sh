#!/usr/bin/env bash
###############################################################################
##  run_all.sh -- runs the empirical example: run_model.r and
##  synthetic_instrument.r for lambda = -0.8, -0.6, -0.4, -0.2, in parallel,
##  inside the pinned Nix environment (see gen-env.R).
##
##  Usage (from the repository root, inside `nix-shell --pure` or via direnv):
##    RESULTS_DIR=emperical_example/results bash emperical_example/run_all.sh
##
##  RESULTS_DIR keeps runs with different rapsem versions apart; every script
##  writes its results and a log file there.
###############################################################################
set -euo pipefail
cd "$(dirname "$0")/.."

export RESULTS_DIR="${RESULTS_DIR:-emperical_example/results}"
export OPENBLAS_NUM_THREADS=1
mkdir -p "$RESULTS_DIR"

run() {
  if [[ "$1" == main ]]; then
    Rscript emperical_example/run_model.r > "$RESULTS_DIR/run_model.log" 2>&1
  else
    LAMBDA="$1" Rscript emperical_example/synthetic_instrument.r \
      > "$RESULTS_DIR/synthetic_instrument_lambda$1.log" 2>&1
  fi
  echo "done $1 ($(date))"
}
export -f run

echo "rapsem $(Rscript -e 'cat(as.character(packageVersion("rapsem")))') -> $RESULTS_DIR ($(date))"
printf '%s\n' main -0.80 -0.60 -0.40 -0.20 | xargs -P 5 -I{} bash -c 'run {} || echo "FAILED {}"'
