#!/usr/bin/env Rscript
# Choose the gradient forest predictors for one species.
# Candidates = the requested layers that vary across the populations with genotype
# data. Among correlated candidates (|Pearson r| > threshold at those populations)
# keep the one ranked higher by rank_by (see select_predictors() in gf_common.R);
# with rank_by importance/accuracy that is the variable that fits the allele
# frequencies better in a screening fit with all candidates. The randomized fits
# use the selected set.
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))
suppressPackageStartupMessages(library(ggplot2))

opt <- parse_args()
need(opt, c("env_sites", "popfreq", "layers", "include", "threshold", "rank_by", "keep_priority",
            "screen_ntree", "maxLevel", "corr_threshold", "nbin", "seed",
            "out_list", "out_candidates", "out_report", "out_heatmap"))
if (!opt$rank_by %in% c("importance", "accuracy", "correlation"))
  stop("rank_by must be importance, accuracy or correlation, got ", opt$rank_by)
threshold <- if (opt$threshold == "none") NA_real_ else as.numeric(opt$threshold)
keep <- if (opt$keep_priority == "none") character(0) else strsplit(opt$keep_priority, ",")[[1]]
floor <- parse_floor(opt, opt$rank_by)
if (!is.na(floor) && is.na(threshold)) stop("importance_floor needs predictor selection (a corr_threshold)")

layers <- trimws(readLines(opt$layers))
layers <- layers[nzchar(layers)]
vars <- if (opt$include == "all") layers else strsplit(opt$include, ",")[[1]]
unknown <- setdiff(c(vars, keep), layers)
if (length(unknown)) stop("predictors not among layers: ", paste(unknown, collapse = ", "))

# load every requested layer for the populations with genotype data
writeLines(vars, opt$out_candidates)
d <- load_gf_data(opt$env_sites, opt$popfreq, opt$out_candidates)
const <- vars[vapply(d$X, function(x) isTRUE(sd(x) == 0), logical(1))]
cand <- setdiff(vars, const)
writeLines(cand, opt$out_candidates)
X <- d$X[, cand, drop = FALSE]

set.seed(as.integer(opt$seed))
sel <- select_predictors(X, d$Y, opt, threshold, opt$rank_by, keep, floor)

score_col <- switch(opt$rank_by, importance = "screen_importance", accuracy = "screen_accuracy",
                    correlation = "screen_score")
rep <- data.table(variable = cand,
                  status = ifelse(cand %in% sel$selected, "kept", "dropped"),
                  rank = sel$rank[cand],
                  score = signif(sel$score[cand], 4),
                  screen_share = round(sel$share[cand], 4),
                  mean_abs_r = round(sel$mean_abs_r[cand], 3),
                  reason = sel$reason[cand])
setnames(rep, "score", score_col)
if (opt$rank_by == "correlation") rep[, c(score_col, "screen_share") := NULL]
setorder(rep, rank)
if (length(const)) rep <- rbind(rep, data.table(variable = const, status = "dropped",
                                                reason = "constant across sites"), fill = TRUE)
fwrite(rep, opt$out_report, sep = "\t", na = "NA", quote = FALSE)
writeLines(sel$selected, opt$out_list)
message(sprintf("%d sites; kept %d of %d predictors (rank_by = %s, threshold = %s, importance_floor = %s)",
                nrow(X), length(sel$selected), length(vars), opt$rank_by, opt$threshold,
                if (is.na(floor)) "none" else floor))
print(rep)

# correlation heatmap of the candidates; kept variables in bold black
cm <- cor(X)
long <- data.table(a = rep(cand, times = length(cand)), b = rep(cand, each = length(cand)), r = as.vector(cm))
long[, a := factor(a, levels = cand)][, b := factor(b, levels = rev(cand))]
is_kept <- cand %in% sel$selected
p <- ggplot(long, aes(a, b, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = "Pearson r") +
  labs(x = NULL, y = NULL,
       title = sprintf("Predictor correlations at %d sites (bold = kept)", nrow(X))) +
  theme_minimal(base_size = 10) +
  # vectorized element_text() works but warns that it is not officially supported
  suppressWarnings(theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   face = ifelse(is_kept, "bold", "plain"),
                                   colour = ifelse(is_kept, "black", "grey50")),
        axis.text.y = element_text(face = rev(ifelse(is_kept, "bold", "plain")),
                                   colour = rev(ifelse(is_kept, "black", "grey50"))),
        panel.grid = element_blank()))
sz <- 0.42 * length(cand) + 2
suppressWarnings(ggsave(opt$out_heatmap, p, width = sz + 1, height = sz, dpi = 150, bg = "white"))
