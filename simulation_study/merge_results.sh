#!/usr/bin/env bash
###############################################################################
##  merge_results.sh -- checks and merges the per-cell results into the
##  per-study files that the plotting scripts read:
##    results/study_1/study_1_results.csv, study_1_fit_indices.csv,
##    results/study_2/study_2_results.csv
##
##  Usage (from anywhere; needs only bash):
##    bash simulation_study/merge_results.sh
##  RESULTS_DIR=<dir> writes the merged files below <dir> instead of
##  simulation_study/results (the cell files are always read from
##  simulation_study/results/cells).
##
##  Every cell of cells.sh is checked: finished (<csv>.done), 100 replicates
##  (5 rows each: med for standard/gestimation/lavaan + two condition
##  numbers), and 100 fit-index rows where fit indices are requested. Cells
##  are appended in the design order of cells.sh. Exit status 1 if any cell
##  is missing or incomplete (the complete cells are merged anyway).
###############################################################################
set -euo pipefail
cd "$(dirname "$0")/.."

CELLS="simulation_study/results/cells"
OUT="${RESULTS_DIR:-simulation_study/results}"
REPS=100
ROWS_PER_REP=5

target() {  # <study> <results|fit>
  case "$1" in
    study_1) echo "$OUT/study_1/study_1_$([[ $2 == fit ]] && echo fit_indices || echo results).csv" ;;
    study_2) echo "$OUT/study_2/study_2_results.csv" ;;
  esac
}

mkdir -p "$OUT/study_1" "$OUT/study_2"
declare -A started=()
missing=() incomplete=()
total=0 complete=0

append() {  # <cell file> <target file>
  if [[ -z "${started[$2]:-}" ]]; then
    head -n 1 "$1" > "$2"
    started[$2]=1
  fi
  tail -n +2 "$1" >> "$2"
}

while read -r study name _; do
  total=$((total + 1))
  csv="$CELLS/${study}__${name}.csv"
  fit="$CELLS/${study}__${name}_fit.csv"
  if [[ ! -f "$csv.done" ]]; then
    missing+=("${study}__${name}")
    continue
  fi
  rows=$(($(wc -l < "$csv") - 1))
  ok=yes
  if (( rows != REPS * ROWS_PER_REP )); then
    incomplete+=("${study}__${name}: $((rows / ROWS_PER_REP)) of $REPS replicates")
    ok=no
  fi
  append "$csv" "$(target "$study" results)"
  if [[ "$study" == study_1* ]]; then
    fit_rows=$(( $(wc -l < "$fit") - 1 ))
    if (( fit_rows != REPS )); then
      incomplete+=("${study}__${name}: $fit_rows of $REPS fit-index rows")
      ok=no
    fi
    append "$fit" "$(target "$study" fit)"
  fi
  [[ "$ok" == yes ]] && complete=$((complete + 1))
done < <(bash simulation_study/cells.sh)

for f in "${!started[@]}"; do
  echo "merged: $f ($(($(wc -l < "$f") - 1)) rows)"
done | sort
echo "cells complete: $complete of $total"
if (( ${#missing[@]} > 0 )); then
  echo "missing or unfinished cells (${#missing[@]}):"
  printf '  %s\n' "${missing[@]}"
fi
if (( ${#incomplete[@]} > 0 )); then
  echo "incomplete cells (see results/cells/logs/ for the failed replicates):"
  printf '  %s\n' "${incomplete[@]}"
fi
(( complete == total ))
