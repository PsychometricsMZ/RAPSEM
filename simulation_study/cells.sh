#!/usr/bin/env bash
###############################################################################
##  cells.sh -- the cells of the simulation study, one per line:
##    <study> <name> <num_obs> <interaction> <theta_m> <rel> <conf> <conf_med> <with_fit>
##
##  Usage:  bash simulation_study/cells.sh [study]
##          study = study_1 | study_2 (all cells if omitted)
##          study 1: confounding (conf_med = 0) and confounded mediator
##          (conf_med > 0); study 2: power over theta, interaction,
##          reliability and n
##
##  Single source of the design for run_local.sh (local, parallel), the
##  run_study_{1,2}.sh SLURM array jobs (one array task per cell) and
##  merge_results.sh. The study label also selects the merged result file.
###############################################################################
set -euo pipefail

all_cells() {
  local n conf conf_med theta inter rel
  for n in 250 500 750 1000; do
    for conf in 0.0 0.100 0.245 0.316 0.374; do
      echo "study_1 n${n}_conf${conf} $n 0.232 0.0 0.9 $conf 0.0 yes"
    done
    for conf_med in 0.100 0.245 0.316 0.374; do
      echo "study_1 n${n}_confmed${conf_med} $n 0.232 0.0 0.9 0.374 $conf_med yes"
    done
  done
  for n in 250 500 750 1000; do
    for theta in 0.0 0.100 0.243 0.316 0.374; do
      for inter in 0.062 0.152 0.196 0.232; do
        for rel in 0.5 0.7 0.8 0.9; do
          echo "study_2 n${n}_theta${theta}_inter${inter}_rel${rel} $n $inter $theta $rel 0.0 0.0 no"
        done
      done
    done
  done
}

if [[ $# -eq 0 ]]; then
  all_cells
else
  all_cells | awk -v s="$1" 'index($1, s) == 1'
fi
