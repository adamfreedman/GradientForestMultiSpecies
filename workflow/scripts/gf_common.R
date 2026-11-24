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
  # gradientForest fits one randomForest per SNP inside try(): collect the
  # warnings and the try() error text instead of letting them scroll past.
  warns <- character(0)
  err_lines <- character(0)
  tc <- textConnection("err_lines", "w", local = TRUE)
  old <- options(try.outFile = tc)
  on.exit({ options(old); close(tc) }, add = TRUE)
  gf <- withCallingHandlers(
    gradientForest(cbind(X, Y),
                   predictor.vars = colnames(X), response.vars = colnames(Y),
                   ntree = as.integer(opt$ntree), transform = NULL, compact = FALSE,
                   nbin = as.integer(opt$nbin), maxLevel = maxLevel,
                   corr.threshold = as.numeric(opt$corr_threshold), trace = FALSE),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  n_err <- sum(grepl("^Error", err_lines))
  attr(gf, "maxLevel") <- maxLevel
  attr(gf, "input") <- list(n_input = ncol(Y), n_pos_rsq = as.integer(gf$species.pos.rsq),
                            n_errors = n_err, warnings = table(warns))
  message(sprintf("gradientForest input: %d SNPs; entered the model (R^2 > 0): %d; randomForest errors: %d",
                  ncol(Y), gf$species.pos.rsq, n_err))
  if (length(warns)) {
    tw <- sort(table(warns), decreasing = TRUE)
    for (i in seq_along(tw)) message(sprintf("  randomForest warning x %d: %s", tw[[i]], names(tw)[i]))
  }
  if (n_err) message("  first error: ", paste(head(err_lines, 2), collapse = " "))
  gf
}

# Input SNPs by what randomForest can do with them.
# few_unique: <= 5 distinct population frequencies (randomForest warns that a
# regression on these may not be meant); constant: one value, so R^2 is NA and
# the SNP can never enter the model.
snp_input_counts <- function(Y) {
  nu <- vapply(Y, function(y) length(unique(y[!is.na(y)])), integer(1))
  list(n_few_unique = sum(nu <= 5), n_constant = sum(nu <= 1))
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

# Choose predictors among correlated environmental variables.
# Variables are ranked, then taken in rank order, skipping any whose |Pearson r|
# at the sites with a variable already kept is above threshold; so no two kept
# variables are correlated above threshold, and of a correlated set the
# higher-ranked one is kept. rank_by:
#   importance / accuracy: R^2-weighted / accuracy importance in a screening
#     gradient forest fitted with all candidates (screen_ntree trees), i.e. the
#     variable that fits the allele frequencies better is kept;
#   correlation: lowest mean |r| with the other candidates first (ignores the SNPs).
# keep_priority variables are ranked first. threshold NA = keep every candidate.
select_predictors <- function(X, Y, opt, threshold, rank_by, keep_priority = character(0)) {
  v <- colnames(X)
  r <- abs(cor(X))
  diag(r) <- 0
  mean_abs_r <- if (length(v) > 1) colSums(r) / (length(v) - 1) else setNames(0, v)
  score <- setNames(rep(NA_real_, length(v)), v)
  if (!is.na(threshold) && length(v) > 1) {
    if (rank_by == "correlation") {
      score <- -mean_abs_r
    } else {
      o <- opt
      o$ntree <- opt$screen_ntree
      message("screening fit with all ", length(v), " candidate predictors")
      gf <- fit_gf(X, Y, o)
      score <- gf_stats(gf, v)[[rank_by]]
      rm(gf)
      invisible(gc())
    }
  }
  ord <- if (all(is.na(score))) seq_along(v) else order(!(v %in% keep_priority), -score, seq_along(v))
  kept <- character(0)
  reason <- setNames(rep("", length(v)), v)
  for (x in v[ord]) {
    if (is.na(threshold) || !length(kept) || max(r[x, kept]) <= threshold) {
      kept <- c(kept, x)
    } else {
      w <- kept[which.max(r[x, kept])]
      reason[x] <- sprintf("|r|=%.3f with %s, which ranks higher", r[x, w], w)
    }
  }
  list(selected = v[v %in% kept], rank = setNames(match(v, v[ord]), v),
       score = if (rank_by == "correlation") setNames(rep(NA_real_, length(v)), v) else score,
       mean_abs_r = mean_abs_r, reason = reason)
}
