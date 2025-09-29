#!/bin/bash
#SBATCH --job-name=study1
#SBATCH --output=logs/output_%j.out
#SBATCH --error=logs/errors_%j.err
#SBATCH --mem=2000
#SBATCH --time=4-00:00:00
#SBATCH --partition=cpu-single
#SBATCH --nodes=1

module purge
module load math/R/4.4.2

cd "$SLURM_SUBMIT_DIR"

output_file="results/study_1/study_1_results.csv"

Rscript init_file.R $output_file

add_var=1e-4
beta_mxr=0.204
beta_ym=0.0
rel=0.8

for num_obs in 100 250 500 750 1000; do
  for conf in 0.0 0.2 0.4 0.6; do
    for conf_treat in 0.0 0.3 0.6 0.9; do
      Rscript run_simulation.R "$add_var" "$num_obs" "$beta_mxr" "$beta_ym" "$rel" "$conf" "$conf_treat" $output_file
    done
  done
done
