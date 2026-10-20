# RAPSEM: Rank-Preserving Structural Equation Modeling

Code accompanying the paper [arXiv:2509.23935](https://arxiv.org/abs/2509.23935)
(R package, simulation studies, and empirical example).

RAPSEM estimates latent mediation effects in randomized controlled trials
that are robust to unobserved mediator–outcome confounding. It combines
factor scores with G-estimation, using treatment × covariate interactions as
instruments, and corrects for measurement error with the method of moments of
Wall and Amemiya (2000). Models are specified in `lavaan` syntax.

## Repository

| Path | Content |
|---|---|
| `rapsem/` | R package (source); built as `rapsem_0.0.0.2000.tar.gz` |
| `simulation_study/` | Simulation studies 1 (confounding) and 2 (power) |
| `emperical_example/` | Empirical example (analysis, tables, figure) |
| `default.nix`, `gen-env.R`, `docker-image.nix` | Pinned R environment (Nix/rix) and container image |

## Installation

```r
install.packages(c("R6", "lavaan"))
install.packages("rapsem_0.0.0.2000.tar.gz", repos = NULL, type = "source")
```

```r
library(rapsem)
model <- "
  cov =~ x1 + x2 + x3
  med =~ m1 + m2 + m3
  out =~ y1 + y2 + y3
  med ~ treat + cov + treat:cov
  out ~ treat + med + cov
"
fit <- med_est(model, data, treat_names = "treat")
fit$summary$gestimation   # RAPSEM (G-estimation)
fit$summary$standard      # factor score regression without G-estimation
```

Estimates are unstandardized by default (metric of the first indicator of
each latent variable); `standardize = TRUE` reports standardized estimates.
Standard errors and percentile confidence intervals are obtained by
bootstrapping (`nboot = 100`).

## Reproducing the results

The analyses were run in the pinned environment of `default.nix`
(R 4.6.1, lavaan 0.7.2; `nix-shell --pure`). Without Nix, the scripts
additionally need the R packages here, mvtnorm, dplyr, tidyr, stringr, haven,
and mice. All commands are run from the repository root.

**Simulation studies**

```bash
bash simulation_study/run_local.sh                # locally, in parallel (JOBS=12)
sbatch simulation_study/run_study_1.sh            # or on a SLURM cluster,
sbatch simulation_study/run_study_2.sh            # with rapsem.sif in simulation_study/
bash simulation_study/merge_results.sh            # merge the per-cell results
```

Figures: `python simulation_study/results/study_1/*.py` and
`python simulation_study/results/study_2/*.py` (pandas, numpy, matplotlib,
statsmodels), written to `simulation_study/results/study_*/plots/`.

**Empirical example**

```bash
bash emperical_example/run_all.sh                 # results in emperical_example/results/
python emperical_example/table.py                 # mediation table
python emperical_example/table_instrument.py      # synthetic instrument table
Rscript emperical_example/table_reliability.r     # reliability table
python emperical_example/plot.py                  # figure
```

## License

MIT (see `rapsem/LICENSE`).
