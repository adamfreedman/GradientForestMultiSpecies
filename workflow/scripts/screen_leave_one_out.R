#!/usr/bin/env Rscript
# Predictor selection with one population left out: a diagnostic of how much
# the choice of predictors depends on which populations were sampled.
# Repeats select_predictors() (screening fit + correlation rule) on the
# populations with genotype data minus population number --k; correlations are
# recomputed at the remaining sites. Writes one row per candidate variable.
# If --k is larger than the number of populations, writes a header-only file.
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))

opt <- parse_args()
need(opt, c("env_sites", "popfreq", "candidates", "threshold", "rank_by", "keep_priority",
            "screen_ntree", "maxLevel", "corr_threshold", "nbin", "seed", "k", "out"))
threshold <- if (opt$threshold == "none") NA_real_ else as.numeric(opt$threshold)
keep <- if (opt$keep_priority == "none") character(0) else strsplit(opt$keep_priority, ",")[[1]]
k <- as.integer(opt$k)

empty <- data.table(left_out = character(0), sitename = character(0), variable = character(0),
                    score = numeric(0), rank = integer(0), kept = logical(0),
                    nearest_kept = character(0), r_nearest = numeric(0), n_pos_rsq = integer(0))
d <- load_gf_data(opt$env_sites, opt$popfreq, opt$candidates)
if (k > nrow(d$X)) {
  message(sprintf("k = %d > %d populations: nothing to leave out", k, nrow(d$X)))
  fwrite(empty, opt$out, sep = "\t")
  quit(status = 0)
}
left_out <- d$sites$popid[k]
X <- d$X[-k, , drop = FALSE]
Y <- d$Y[-k, , drop = FALSE]
# a variable can become constant without this site
const <- names(X)[vapply(X, function(x) isTRUE(sd(x) == 0), logical(1))]
X <- X[, setdiff(names(X), const), drop = FALSE]
message(sprintf("leaving out population %s (%s): %d sites, %d candidates%s", left_out, d$sites$sitename[k],
                nrow(X), ncol(X), if (length(const)) paste0("; constant without it: ", paste(const, collapse = ", ")) else ""))

set.seed(as.integer(opt$seed) + k)
sel <- select_predictors(X, Y, opt, threshold, opt$rank_by, intersect(keep, names(X)))

# for every variable, the most correlated kept variable (other than itself):
# for a variable dropped here, the one that stood in for it
r <- abs(cor(X))
v <- names(X)
others <- lapply(v, function(x) setdiff(sel$selected, x))
nearest <- vapply(seq_along(v), function(i) if (length(others[[i]])) others[[i]][which.max(r[v[i], others[[i]]])] else NA_character_, "")
r_near <- vapply(seq_along(v), function(i) if (is.na(nearest[i])) NA_real_ else r[v[i], nearest[i]], 0)

out <- data.table(left_out = left_out, sitename = d$sites$sitename[k], variable = v,
                  score = sel$score[v], rank = sel$rank[v], kept = v %in% sel$selected,
                  nearest_kept = nearest, r_nearest = round(r_near, 3), n_pos_rsq = sel$n_pos_rsq)
if (length(const)) out <- rbind(out, data.table(left_out = left_out, sitename = d$sites$sitename[k],
                                                variable = const, kept = FALSE, n_pos_rsq = sel$n_pos_rsq), fill = TRUE)
fwrite(out, opt$out, sep = "\t", na = "NA")
message("kept: ", paste(sel$selected, collapse = ", "))
