###############################################################################
##  gen-env.R -- generates default.nix: R, Python, and every package version
##  for this project (rapsem + simulation_study + emperical_example), pinned
##  via Nix/rix so the environment stays identical across machines and time.
##
##  One-time setup (needs WSL2 on Windows -- see
##  https://docs.ropensci.org/rix/articles/setting-up-linux-windows.html):
##    curl --proto '=https' --tlsv1.2 -sSf -L https://install.determinate.systems/nix | sh -s -- install
##
##  Regenerate default.nix (only when the environment itself changes), from
##  the repo root inside WSL:
##    nix-shell -I nixpkgs=https://github.com/rstats-on-nix/nixpkgs/archive/2026-08-24.tar.gz \
##      -p R rPackages.rix --run "Rscript gen-env.R"
##
##  Build and enter the environment:
##    nix-build          # builds the environment (first run downloads everything)
##    nix-shell --pure   # enters it -- `type R` should print a /nix/store/... path
##
##  direnv (.envrc, `direnv allow` once) activates this automatically on cd,
##  including in an editor's integrated terminal.
##
##  Inside the shell, `rapsem` is already installed from the bundled
##  rapsem_0.0.0.2000.tar.gz. After changing anything under rapsem/R/:
##  rebuild the tarball (R CMD build rapsem), update local_r_pkgs below if the
##  version changed, re-run this file and nix-build.
##
##  For running on an HPC cluster via Apptainer, see docker-image.nix.
###############################################################################

library(rix)

project_directory <- "."

rix(
  date = "2026-08-24",   # pins R + every CRAN package to this rstats-on-nix snapshot
  r_pkgs = c(
    "R6", "lavaan",      # rapsem/DESCRIPTION: Imports (MASS ships with R)
    "here", "mvtnorm",   # simulation_study/data_simulation
    "dplyr", "haven", "stringr", "mice", "tidyr",  # emperical_example
    "devtools",          # build/test rapsem; builds the .tar.gz below
    "roxygen2",          # regenerate NAMESPACE/man/ (rapsem/DESCRIPTION: RoxygenNote)
    "testthat"           # rapsem's test suite
  ),
  local_r_pkgs = "rapsem_0.0.0.2000.tar.gz",  # built via R CMD build rapsem
  py_conf = list(
    py_version = "3.12",
    py_pkgs    = c("pandas", "numpy", "matplotlib", "statsmodels")  # plots and tables
  ),
  system_pkgs = "apptainer",  # container runtime -- matches many HPC clusters' setup
  ide = "none",   # keep using the system's own VS Code, via direnv (.envrc)
  project_path = project_directory,
  overwrite = TRUE,
  print = TRUE
)
