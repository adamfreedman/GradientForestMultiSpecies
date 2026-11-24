#!/usr/bin/env Rscript
# Randomization replicates rep_start..rep_end: permute whole rows of the
# environment table (so predictor correlations are kept but decoupled from
# allele frequencies) and refit with the same settings as the observed model.
# Each replicate has its own seed, so results do not depend on batch size.
# The predictors are the observed model's own set, so each SNP's observed R^2 is
# compared with what permuted environments give under the same model.
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))

opt <- parse_args()
need(opt, c("env_sites", "popfreq", "predictors", "ntree_file", "maxLevel", "corr_threshold", "nbin",
            "seed", "rep_start", "rep_end", "out", "out_snp_rsq"))
opt$ntree <- fread(opt$ntree_file)$ntree  # same number of trees as the observed model
d <- load_gf_data(opt$env_sites, opt$popfreq, opt$predictors)

one_rep <- function(rep) {
  set.seed(as.integer(opt$seed) * 100003L + rep)
  perm <- sample(nrow(d$X))
  Xr <- d$X[perm, , drop = FALSE]
  rownames(Xr) <- NULL
  gf <- fit_gf(Xr, d$Y, opt)
  st <- gf_stats(gf, d$preds)
  # per-SNP R^2 (only SNPs with R^2 > 0 are reported by gradientForest; absent = R^2 <= 0)
  snp_rsq[[length(snp_rsq) + 1]] <<- if (is.null(gf)) NULL else
    data.table(rep = rep, snp = names(gf$result), rsq = signif(unname(gf$result), 5))
  row <- data.table(rep = rep, perm = paste(perm, collapse = ","),
                    n_pos_rsq = st$n_pos_rsq, mean_rsq = st$mean_rsq, sum_rsq = st$sum_rsq)
  cbind(row,
        as.data.table(as.list(setNames(st$importance, paste0("imp_", names(st$importance))))),
        as.data.table(as.list(setNames(st$accuracy, paste0("acc_", names(st$accuracy))))))
}

snp_rsq <- list()
reps <- seq(as.integer(opt$rep_start), as.integer(opt$rep_end))
res <- rbindlist(lapply(reps, function(r) {
  t0 <- Sys.time()
  x <- one_rep(r)
  message(sprintf("rep %d done in %.1f min", r, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  x
}))
fwrite(res, opt$out, sep = "\t")
fwrite(rbindlist(snp_rsq), opt$out_snp_rsq, sep = "\t")  # .gz suffix -> gzip
print(res[, !"perm"])
