#!/bin/bash
#SBATCH --job-name=study2
#SBATCH --output=logs/output_%j.out
#SBATCH --error=logs/errors_%j.err
#SBATCH --mem=2000
#SBATCH --time=4-00:00:00
#SBATCH --partition=cpu-single
#SBATCH --nodes=1

# Load R module
module purge
module load math/R/4.4.2

cd "$SLURM_SUBMIT_DIR"

output_file="results/study_2/study_2_results.csv"

add_var=1e-4
conf=0.4
conf_treat=0.0

for num_obs in 1000; do
  for beta_mxr in 0.102 0.145 0.176 0.204; do
    for beta_ym in 0.41; do
      for rel in 0.667 0.8; do
        Rscript run_simulation.R "$add_var" "$num_obs" "$beta_mxr" "$beta_ym" "$rel" "$conf" "$conf_treat" $output_file
      done
    done
  done
done

for num_obs in 1000; do
  for beta_mxr in 0.145 0.176 0.204; do
    for beta_ym in 0.0 0.29 0.41; do
      for rel in 0.4 0.667 0.8; do
        Rscript run_simulation.R "$add_var" "$num_obs" "$beta_mxr" "$beta_ym" "$rel" "$conf" "$conf_treat" $output_file
      done
    done
  done
done
