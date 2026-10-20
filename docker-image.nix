###############################################################################
##  docker-image.nix -- packages the same pinned environment as default.nix
##  (see gen-env.R for how that's built/entered directly) into an OCI image,
##  for running on an HPC cluster via Apptainer (many HPC clusters support
##  only Apptainer/Singularity, not Nix or Docker directly).
##
##  Build:
##    nix-build docker-image.nix -o result-image
##
##  Convert to a Singularity image, no local Docker needed. dockerTools'
##  output is gzip-compressed as a whole (docker load auto-detects that) but
##  Apptainer's docker-archive reader expects an uncompressed outer tar
##  (compression belongs to individual layers, not the envelope) -- gunzip
##  first or the build fails with "gzip: invalid header":
##    gunzip -k result-image -c > image.tar   # or: gunzip -k result-image (renames in place)
##    apptainer build rapsem.sif docker-archive://image.tar
##
##  On the cluster (or apptainer exec, if that's what's available there):
##    singularity exec rapsem.sif Rscript simulation_study/run_simulation.R ...
##
##  Rebuild whenever default.nix changes (new gen-env.R run + nix-build).
###############################################################################

let
  env = import ./default.nix;
  pkgs = env.pkgs;
  allInputs = pkgs.lib.flatten env.shell.buildInputs;

  # nix-shell wires up R_LIBS_SITE via each R package's setup-hook at
  # shell-entry time, walking the *transitive* dependency graph (fs, cli,
  # rlang, ... -- not just the directly-requested packages). A plain
  # container runtime never runs that hook, and reimplementing the closure
  # walk by hand only recovers the direct list, not transitive deps. Instead,
  # let rWrapper do what it already does: build a single wrapped R/Rscript
  # whose transitive R_LIBS_SITE is resolved and baked in at build time.
  isRRelated = p: let n = p.name or ""; in
    pkgs.lib.hasPrefix "R-" n || pkgs.lib.hasPrefix "r-" n;
  rPkgs = pkgs.lib.filter isRRelated allInputs;
  rEnv = pkgs.rWrapper.override { packages = rPkgs; };

  # nixpkgs python packages need the same treatment: wired into a single
  # interpreter via withPackages (site-packages layout), not just present in
  # the store -- a bare `python3` plus separate pandas/matplotlib derivations
  # on PATH doesn't make them importable. Distinguish site-packages
  # ("python3.12-pandas-3.0.4", note the dot) from the bare interpreter
  # ("python3-3.12.14", note the dash) by name prefix.
  isPythonPkg = p: pkgs.lib.hasPrefix "python3.12-" (p.name or "");
  isBarePython = p: pkgs.lib.hasPrefix "python3-" (p.name or "");
  pyPkgs = pkgs.lib.filter isPythonPkg allInputs;
  pyEnv = pkgs.python312.withPackages (ps: pyPkgs);

  # apptainer and nix itself are in default.nix's system_pkgs purely for
  # local dev convenience (building/testing this very image) -- pointless
  # and bloating to ship inside the image that Apptainer runs on the cluster.
  isLocalOnlyTool = p: builtins.elem (p.pname or "") [ "apptainer" "nix" ];
  otherInputs = pkgs.lib.filter
    (p: !isRRelated p && !isLocalOnlyTool p && !isPythonPkg p && !isBarePython p)
    allInputs;
in
pkgs.dockerTools.buildImage {
  name = "rapsem-simulation-study";
  tag = "latest";

  copyToRoot = pkgs.buildEnv {
    name = "image-root";
    paths = otherInputs ++ [
      rEnv
      pyEnv
      pkgs.bashInteractive
      pkgs.coreutils
    ];
    pathsToLink = [ "/bin" "/lib" "/share" ];
  };

  config = {
    Env = [
      "LANG=en_US.UTF-8"
      "LC_ALL=en_US.UTF-8"
      "LOCALE_ARCHIVE=${pkgs.glibcLocales}/lib/locale/locale-archive"
    ];
    Cmd = [ "${pkgs.bashInteractive}/bin/bash" ];
    WorkingDir = "/workspace";
  };
}
