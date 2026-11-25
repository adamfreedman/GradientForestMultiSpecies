#!/usr/bin/env Rscript
# Fit the observed gradient forest (population allele frequencies ~ environment).
#
# With --ntree=auto the number of trees is chosen adaptively: fit at
# tune_start trees, then multiply by tune_factor and refit (new seed each time)
# until two successive fits agree on
#   * per-SNP R^2            (Pearson r >= min_snp_rsq_cor; SNPs absent from $result count as 0)
#   * predictor importance   (max |change| in R^2-weighted importance <= max_importance_change,
#                             as a fraction of total importance)
#   * turnover space         (Pearson r of pairwise distances between n_cells random grid
#                             cells in transformed-predictor space >= min_turnover_cor)
# or tune_max trees is reached (then a warning is written to the log and to
# ntree_selected.tsv). The larger fit of the first agreeing pair is kept.
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))
suppressPackageStartupMessages(library(terra))

opt <- parse_args()
need(opt, c("env_sites", "popfreq", "predictors", "env_tif", "ntree", "maxLevel", "corr_threshold",
            "nbin", "seed", "out_model", "out_importance", "out_rsq", "out_summary",
            "out_tuning", "out_ntree"))
seed <- as.integer(opt$seed)
d <- load_gf_data(opt$env_sites, opt$popfreq, opt$predictors)
snps <- colnames(d$Y)

fit_at <- function(ntree, k) {
  set.seed(seed + k - 1L)
  o <- opt
  o$ntree <- ntree
  t0 <- Sys.time()
  gf <- fit_gf(d$X, d$Y, o)
  if (is.null(gf)) stop("gradientForest returned NULL at ntree=", ntree, ": no SNP had positive R^2")
  attr(gf, "minutes") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  gf
}

auto <- identical(opt$ntree, "auto")
tuning <- data.table()
warn_msg <- ""
if (!auto) {
  ntree <- as.integer(opt$ntree)
  gf <- fit_at(ntree, 1L)
  converged <- NA
  tuning <- data.table(step = 1L, ntree = ntree, minutes = attr(gf, "minutes"),
                       n_pos_rsq = gf$species.pos.rsq)
} else {
  need(opt, c("tune_start", "tune_max", "tune_factor", "min_snp_rsq_cor",
              "max_importance_change", "min_turnover_cor", "n_cells"))
  t_start <- as.integer(opt$tune_start)
  t_max <- as.integer(opt$tune_max)
  t_fac <- as.numeric(opt$tune_factor)
  thr_rsq <- as.numeric(opt$min_snp_rsq_cor)
  thr_imp <- as.numeric(opt$max_importance_change)
  thr_turn <- as.numeric(opt$min_turnover_cor)
  if (t_fac <= 1) stop("tune_factor must be > 1")
  if (t_start > t_max) stop("tune_start must be <= tune_max")

  # fixed random sample of grid cells for comparing turnover space
  env <- rast(opt$env_tif)[[d$preds]]
  cells <- spatSample(env, size = as.integer(opt$n_cells), method = "random", na.rm = TRUE,
                      as.df = TRUE, values = TRUE)
  cells <- as.data.frame(cells)[, d$preds, drop = FALSE]
  message(nrow(cells), " grid cells sampled for turnover comparison")

  rsq_vec <- function(gf) { v <- setNames(numeric(length(snps)), snps); v[names(gf$result)] <- gf$result; v }
  imp_vec <- function(gf) gf_stats(gf, d$preds)$importance
  turn_dist <- function(gf) as.vector(dist(as.matrix(predict(gf, cells)[, d$preds])))

  schedule <- t_start
  while (tail(schedule, 1) < t_max) schedule <- c(schedule, min(t_max, ceiling(tail(schedule, 1) * t_fac)))

  prev <- NULL
  converged <- FALSE
  for (k in seq_along(schedule)) {
    nt <- schedule[k]
    gf <- fit_at(nt, k)
    cur <- list(rsq = rsq_vec(gf), imp = imp_vec(gf), dist = turn_dist(gf))
    row <- data.table(step = k, ntree = nt, minutes = attr(gf, "minutes"),
                      n_pos_rsq = gf$species.pos.rsq, mean_rsq = mean(gf$result),
                      snp_rsq_cor = NA_real_, importance_max_change = NA_real_,
                      turnover_dist_cor = NA_real_, agrees_with_previous = NA)
    if (!is.null(prev)) {
      row$snp_rsq_cor <- cor(cur$rsq, prev$rsq)
      row$importance_max_change <- max(abs(cur$imp - prev$imp)) / sum(prev$imp)
      row$turnover_dist_cor <- cor(cur$dist, prev$dist)
      row$agrees_with_previous <- row$snp_rsq_cor >= thr_rsq &&
        row$importance_max_change <= thr_imp && row$turnover_dist_cor >= thr_turn
    }
    tuning <- rbind(tuning, row)
    message(sprintf("ntree=%d (%.1f min): snp_rsq_cor=%.4f importance_max_change=%.4f turnover_dist_cor=%.4f",
                    nt, row$minutes, row$snp_rsq_cor, row$importance_max_change, row$turnover_dist_cor))
    if (isTRUE(row$agrees_with_previous)) { converged <- TRUE; break }
    prev <- cur
  }
  ntree <- tail(tuning$ntree, 1)
  if (!converged) {
    warn_msg <- sprintf(paste0("ntree did not converge: results at ntree=%d still differ from ntree=%d ",
                               "(thresholds: snp_rsq_cor >= %s, importance_max_change <= %s, ",
                               "turnover_dist_cor >= %s). Using the maximum, ntree=%d."),
                        ntree, tuning$ntree[max(1, nrow(tuning) - 1)], thr_rsq, thr_imp, thr_turn, ntree)
    warning(warn_msg, call. = FALSE)
    message("WARNING: ", warn_msg)
  }
}
fwrite(tuning, opt$out_tuning, sep = "\t")
fwrite(data.table(ntree = ntree, mode = if (auto) "auto" else "fixed", converged = converged,
                  tune_max = if (auto) as.integer(opt$tune_max) else NA_integer_,
                  warning = warn_msg), opt$out_ntree, sep = "\t")

attr(gf, "popid") <- d$sites$popid  # row order of gf$X
saveRDS(gf, opt$out_model)

st <- gf_stats(gf, d$preds)
imp <- data.table(variable = names(st$importance), importance = st$importance,
                  accuracy = st$accuracy[names(st$importance)])
setorder(imp, -importance)
fwrite(imp, opt$out_importance, sep = "\t")

fwrite(data.table(snp = names(gf$result), rsq = gf$result)[order(-rsq)], opt$out_rsq, sep = "\t")

inp <- attr(gf, "input")
cnt <- snp_input_counts(d$Y)
summ <- data.table(n_sites = nrow(d$X), n_predictors = ncol(d$X), n_snps = ncol(d$Y),
                   n_snps_few_unique = cnt$n_few_unique, n_snps_constant = cnt$n_constant,
                   n_rf_errors = inp$n_errors, n_rf_warnings = sum(inp$warnings),
                   n_pos_rsq = st$n_pos_rsq, mean_rsq = st$mean_rsq, sum_rsq = st$sum_rsq,
                   ntree = ntree, ntree_mode = if (auto) "auto" else "fixed", ntree_converged = converged,
                   maxLevel = attr(gf, "maxLevel"), corr_threshold = as.numeric(opt$corr_threshold))
fwrite(summ, opt$out_summary, sep = "\t")
print(tuning)
print(summ)
print(imp)
