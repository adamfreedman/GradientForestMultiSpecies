# Install extendedForest + gradientForest (R-Forge only; not on CRAN or conda)
# into the first library of the active R (the Snakemake conda env).
lib <- .libPaths()[1]
pkgs <- c("extendedForest", "gradientForest")
have <- function() all(pkgs %in% rownames(installed.packages(lib.loc = lib)))

if (!have()) {
  message("installing from R-Forge into ", lib)
  try(install.packages(pkgs, lib = lib, type = "source",
                       repos = c("https://R-Forge.R-project.org", "https://cloud.r-project.org")))
}
if (!have()) {
  message("R-Forge install failed; trying the GitHub mirror of the R-Forge SVN")
  for (p in pkgs) {
    remotes::install_github("rforge/gradientforest", subdir = file.path("pkg", p),
                            lib = lib, upgrade = "never")
  }
}
suppressPackageStartupMessages(library(gradientForest))
cat("extendedForest", as.character(packageVersion("extendedForest")),
    "gradientForest", as.character(packageVersion("gradientForest")), "\n")
