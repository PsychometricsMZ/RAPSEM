#!/bin/bash
###############################################################################
##  run_study_1.sh -- study 1 of the simulation study (confounding and confounded mediator; theta = 0) as a
##  SLURM array job, one array task per cell (36 cells, 100 replicates each),
##  run inside the pinned container rapsem.sif (R 4.6.1, rapsem 0.0.0.2000;
##  see docker-image.nix for how it is built).
##
##  Usage (from the repository root, with rapsem.sif copied to simulation_study/):
##    sbatch simulation_study/run_study_1.sh
##  When the job has finished, check and merge the per-cell files into the
##  result files read by the plotting scripts (fails if a cell is missing):
##    bash simulation_study/merge_results.sh
##
##  Cells are listed by cells.sh and run by run_cell.sh (both shared with
##  run_local.sh); finished cells are skipped, so failed tasks can simply be
##  resubmitted, e.g. sbatch --array=5,17 simulation_study/run_study_1.sh
###############################################################################
#SBATCH --job-name=rapsem_study1
#SBATCH --output=simulation_study/slurm_logs/study1_%A_%a.out
#SBATCH --error=simulation_study/slurm_logs/study1_%A_%a.err
#SBATCH --array=0-35
#SBATCH --mem=3000
#SBATCH --time=24:00:00
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

set -euo pipefail

module purge
module load system/singularity/3.11.3

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

# rapsem.sif is a build artifact (see docker-image.nix), not committed to git
# -- build it elsewhere (nix-build docker-image.nix; apptainer build) and copy
# it to simulation_study/ before submitting. SingularityCE reads
# Apptainer-built .sif files fine -- same underlying format.
SIF="$(pwd)/simulation_study/rapsem.sif"
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
# --cleanenv: no host environment (e.g. R_LIBS_USER) leaks into the container
export RSCRIPT="singularity exec --cleanenv --bind $(pwd) $SIF Rscript"

cell=$(bash simulation_study/cells.sh study_1 | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
echo "Task ${SLURM_ARRAY_TASK_ID} (job ${SLURM_ARRAY_JOB_ID}): $cell -- $(date)"
# shellcheck disable=SC2086  # the cell line is split into its fields on purpose
bash simulation_study/run_cell.sh $cell
echo "Finished at $(date)"
