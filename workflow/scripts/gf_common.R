# Shared helpers for the gradient forest R scripts.
suppressPackageStartupMessages({
  library(data.table)
})

# --key=value command-line arguments -> named list; bare arguments -> $positional
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(positional = character(0))
  for (a in args) {
    if (grepl("^--[^=]+=", a)) {
      k <- sub("^--([^=]+)=.*$", "\\1", a)
      out[[k]] <- sub("^--[^=]+=", "", a)
    } else {
      out$positional <- c(out$positional, a)
    }
  }
  out
}

need <- function(opt, keys) {
  miss <- setdiff(keys, names(opt))
  if (length(miss)) stop("missing arguments: ", paste0("--", miss, collapse = ", "))
}

# Environment (predictors) and allele frequencies (responses) aligned by popid.
load_gf_data <- function(env_sites, popfreq, predictors) {
  env <- fread(env_sites, colClasses = list(character = "popid"))
  freq <- fread(popfreq, colClasses = list(character = "popid"))
  preds <- readLines(predictors)
  preds <- preds[nzchar(preds)]
  common <- intersect(env$popid, freq$popid)
  dropped <- setdiff(env$popid, common)
  if (length(dropped)) message("populations without genotype data (dropped): ", paste(dropped, collapse = ", "))
  env <- env[match(common, popid)]
  freq <- freq[match(common, popid)]
  X <- as.data.frame(env[, ..preds])
  Y <- as.data.frame(freq[, !"popid"])
  list(X = X, Y = Y, preds = preds, sites = env)
}

resolve_maxlevel <- function(x, n_sites) {
  if (identical(x, "auto")) max(1L, floor(log2(0.368 * n_sites / 2))) else as.integer(x)
}

fit_gf <- function(X, Y, opt) {
  suppressPackageStartupMessages(library(gradientForest))
  maxLevel <- resolve_maxlevel(opt$maxLevel, nrow(X))
  message(sprintf("gradientForest: %d sites, %d predictors, %d SNPs, ntree=%s, maxLevel=%d, corr.threshold=%s",
                  nrow(X), ncol(X), ncol(Y), opt$ntree, maxLevel, opt$corr_threshold))
  gf <- gradientForest(cbind(X, Y),
                       predictor.vars = colnames(X), response.vars = colnames(Y),
                       ntree = as.integer(opt$ntree), transform = NULL, compact = FALSE,
                       nbin = as.integer(opt$nbin), maxLevel = maxLevel,
                       corr.threshold = as.numeric(opt$corr_threshold), trace = FALSE)
  attr(gf, "maxLevel") <- maxLevel
  gf
}

# Summary statistics compared between the real and randomized fits.
# importance = R^2-weighted importance; accuracy = permutation (accuracy) importance,
# averaged over SNPs. Predictors never used get 0.
gf_stats <- function(gf, preds) {
  imp <- acc <- setNames(rep(0, length(preds)), preds)
  if (is.null(gf)) {
    return(list(n_pos_rsq = 0L, mean_rsq = NA_real_, sum_rsq = 0, importance = imp, accuracy = acc))
  }
  w <- importance(gf, type = "Weighted")
  a <- importance(gf, type = "Accuracy")
  imp[names(w)] <- w
  acc[names(a)] <- a
  list(n_pos_rsq = as.integer(gf$species.pos.rsq),
       mean_rsq = mean(gf$result),
       sum_rsq = sum(gf$result),
       importance = imp, accuracy = acc)
}
