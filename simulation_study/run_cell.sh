#!/usr/bin/env bash
###############################################################################
##  run_cell.sh -- runs one cell of the simulation study (100 replicates of
##  run_simulation.R) and writes its own CSV under simulation_study/results/
##  cells/, so parallel cells never append to the same file. Finished cells
##  (marker file <csv>.done) are skipped.
##
##  Usage (from anywhere; arguments = one line of cells.sh):
##    bash simulation_study/run_cell.sh <study> <name> <num_obs> <interaction> \
##      <theta_m> <rel> <conf> <conf_med> <with_fit>
##
##  RSCRIPT selects how R is called, e.g. inside the container on the cluster:
##    RSCRIPT="singularity exec --cleanenv simulation_study/rapsem.sif Rscript"
###############################################################################
set -euo pipefail
cd "$(dirname "$0")/.."

read -r -a rscript <<< "${RSCRIPT:-Rscript}"
study="$1" name="$2" n="$3" inter="$4" theta="$5" rel="$6" conf="$7" conf_med="$8" fit="$9"

cells="simulation_study/results/cells"
mkdir -p "$cells/logs"
out="$cells/${study}__${name}.csv"
fit_out=""
[[ "$fit" == yes ]] && fit_out="$cells/${study}__${name}_fit.csv"

if [[ -f "$out.done" ]]; then
  echo "already done: $out"
  exit 0
fi

"${rscript[@]}" simulation_study/init_file.R "$out" > /dev/null
if [[ -n "$fit_out" ]]; then
  "${rscript[@]}" simulation_study/init_file_model_fit.R "$fit_out" > /dev/null
fi
"${rscript[@]}" simulation_study/run_simulation.R 1e-4 "$n" "$inter" "$theta" "$rel" \
  "$conf" "$out" "$conf_med" "$fit_out" > "$cells/logs/${study}__${name}.log" 2>&1
touch "$out.done"
echo "finished: $out"
