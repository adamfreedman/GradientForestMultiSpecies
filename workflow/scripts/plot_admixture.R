#!/usr/bin/env Rscript
# ADMIXTURE plots (inputs from admixture_summary.py, all in --admixdir):
#   admixture_selected_k  ancestry bar plot at the K with the lowest CV error; when that
#                         K is not stable (pong), K - 1 is shown below it, and when
#                         K - 1 is not stable either, also the largest stable K below
#   admixture_all_k       one panel per K
#   selected_k_modes      at the selected K, the representative run of each pong mode
#   cv_error              CV error by K (every run, and the mean), stable vs not stable
# Individuals are grouped by population (labels as in the PCA plots); populations with
# similar ancestry are placed next to each other, and individuals within a
# population are sorted by its main cluster.
# Cluster colours match across K (clusters are aligned by admixture_summary.py).
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))
suppressPackageStartupMessages(library(ggplot2))

opt <- parse_args()
need(opt, c("species", "admixdir"))
d <- opt$admixdir

theme_set(theme_bw(base_size = 11) +  # black panel border, no gridlines
            theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA)))
# Tableau 10, then its lighter companions: up to 20 clusters
CLUSTER_COLS <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
                  "#9C755F", "#BAB0AC", "#A0CBE8", "#FFBE7D", "#FF9D9A", "#86BCB6", "#8CD17D", "#F1CE63",
                  "#D4A6C8", "#FABFD2", "#D7B5A6", "#79706E")

ksum <- fread(file.path(d, "k_summary.tsv"))
sel <- fread(file.path(d, "selected_k.tsv"))
cv <- fread(file.path(d, "cv_error.tsv"))
modes <- fread(file.path(d, "selected_k_modes.tsv"), colClasses = list(character = c("sample_id", "popid")))
kcv <- sel$K_cv
ks <- ksum$K
read_anc <- function(k) fread(file.path(d, sprintf("ancestry_K%d.tsv", k)),
                              colClasses = list(character = c("sample_id", "popid")))
anc <- setNames(lapply(ks, read_anc), ks)

# individual order: populations with similar ancestry next to each other, then within
# a population by its main cluster. Populations are clustered (average linkage, Euclidean
# distance) on their mean ancestry proportions at every K from 2 to the selected K
# together, so populations that separate only at the highest K stay next to each other;
# the dendrogram's branches are ordered by main cluster at the selected K (then purity).
# The order is the same in every panel and every plot.
base <- copy(anc[[as.character(kcv)]])
qcols <- grep("^cluster", names(base), value = TRUE)
pm <- base[, lapply(.SD, mean), by = popid, .SDcols = qcols]
pop_levels <- pm$popid
num <- suppressWarnings(as.numeric(pop_levels))
pop_levels <- if (all(!is.na(num))) pop_levels[order(num)] else sort(pop_levels)  # PCA plot order
pm <- pm[match(pop_levels, popid)]
M <- as.matrix(pm[, ..qcols])
main_idx <- max.col(M, ties.method = "first")
pop_means <- function(k) {
  a <- anc[[as.character(k)]]
  qc <- grep("^cluster", names(a), value = TRUE)
  m <- a[, lapply(.SD, mean), by = popid, .SDcols = qc]
  as.matrix(m[match(pop_levels, popid), ..qc])
}
if (nrow(M) > 2 && ncol(M) > 1) {
  Mall <- do.call(cbind, lapply(intersect(2:kcv, ks), pop_means))
  dend <- reorder(as.dendrogram(hclust(dist(Mall), method = "average")),
                  wts = main_idx - 0.5 * apply(M, 1, max), agglo.FUN = mean)
  pop_levels <- pop_levels[order.dendrogram(dend)]
}
base[, pop_rank := match(popid, pop_levels)]
base[, main := {
  m <- colMeans(as.matrix(.SD))
  .SD[[which.max(m)]]
}, by = popid, .SDcols = qcols]
setorder(base, pop_rank, -main)
ind_order <- base$sample_id
n_ind <- length(ind_order)
blocks <- base[, .(start = .I[1], end = .I[.N]), by = popid]
blocks[, mid := (start + end) / 2]

long_q <- function(dt, panel) {
  qc <- grep("^cluster", names(dt), value = TRUE)
  x <- melt(dt[, c("sample_id", qc), with = FALSE], id.vars = "sample_id", variable.name = "cluster",
            value.name = "q")
  x[, ind := match(sample_id, ind_order)]
  x[, cluster := factor(cluster, levels = paste0("cluster", seq_along(CLUSTER_COLS)))]
  x[, panel := panel]
  x
}

bar_plot <- function(x, title, subtitle = NULL) {
  x[, panel := factor(panel, levels = unique(panel))]
  rot <- n_ind > 0 && (nrow(blocks) > 12 || max(nchar(blocks$popid)) > 4)
  ggplot(x, aes(ind, q, fill = cluster)) +
    geom_col(width = 1, colour = NA) +
    geom_vline(xintercept = blocks$end[-nrow(blocks)] + 0.5, colour = "white", linewidth = 0.6) +
    facet_wrap(~panel, ncol = 1, strip.position = "top") +
    scale_fill_manual(values = setNames(CLUSTER_COLS, paste0("cluster", seq_along(CLUSTER_COLS))),
                      labels = function(l) sub("cluster", "", l), name = "Cluster", drop = TRUE) +
    scale_x_continuous(breaks = blocks$mid, labels = blocks$popid, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    labs(x = "Population", y = "Ancestry proportion", title = title, subtitle = subtitle) +
    theme(axis.text.x = if (rot) element_text(angle = 90, hjust = 1, vjust = 0.5) else element_text(),
          axis.ticks.x = element_blank(), strip.background = element_rect(fill = "grey92", colour = "black"),
          strip.text = element_text(hjust = 0), panel.spacing = unit(0.4, "lines"),
          legend.key.size = unit(0.8, "lines"))
}
plot_width <- max(7, min(16, 3.5 + 0.04 * n_ind))
save_both <- function(p, name, h) {
  ggsave(file.path(d, paste0(name, ".pdf")), p, width = plot_width, height = h, limitsize = FALSE)
  ggsave(file.path(d, paste0(name, ".png")), p, width = plot_width, height = h, dpi = 300, limitsize = FALSE)
}
k_info <- function(k) {
  r <- ksum[K == k]
  sprintf("%d of %d runs in pong's major mode%s", r$major_mode_runs, r$n_runs,
          if (r$stable == "yes") ", stable" else sprintf(" (%d modes), NOT stable", r$n_modes))
}

# selected K; when it is not stable, also K - 1 and, if K - 1 is not stable either,
# the largest stable K below the selected K
unstable <- sel$stable != "yes"
x <- long_q(anc[[as.character(kcv)]], sprintf("K = %d (lowest CV error): %s", kcv, k_info(kcv)))
n_panels <- 1
if (unstable && (kcv - 1) %in% ks) {
  x <- rbind(x, long_q(anc[[as.character(kcv - 1)]], sprintf("K = %d: %s", kcv - 1, k_info(kcv - 1))))
  n_panels <- 2
}
k_stable <- suppressWarnings(as.integer(sel$largest_stable_K_le_K_cv))
if (unstable && !is.na(k_stable) && k_stable < kcv - 1) {
  x <- rbind(x, long_q(anc[[as.character(k_stable)]],
                       sprintf("K = %d (largest stable K below %d): %s", k_stable, kcv, k_info(k_stable))))
  n_panels <- n_panels + 1
}
save_both(bar_plot(x, sprintf("%s: ADMIXTURE", opt$species),
                   sprintf("%d individuals; each bar is an individual, grouped by population", n_ind)),
          "admixture_selected_k", 1.6 + 2 * n_panels)

# every K
x <- rbindlist(lapply(ks, function(k) long_q(anc[[as.character(k)]], sprintf("K = %d%s%s", k,
  if (k == kcv) " (lowest CV error)" else "", if (ksum[K == k, stable] != "yes") ": not stable" else ""))))
save_both(bar_plot(x, sprintf("%s: ADMIXTURE, K = %d to %d", opt$species, min(ks), max(ks)),
                   "pong major-mode run at each K; cluster colours matched across K"),
          "admixture_all_k", 1.6 + 1.3 * length(ks))

# pong modes at the selected K
MAX_MODES <- 6  # largest modes drawn; all are in selected_k_modes.tsv
mk_all <- unique(modes[, .(mode, rep_run, n_runs, avg_sim_within)])
mk <- head(mk_all, MAX_MODES)
x <- rbindlist(lapply(seq_len(nrow(mk)), function(i) {
  long_q(modes[mode == mk$mode[i], c("sample_id", grep("^cluster", names(modes), value = TRUE)), with = FALSE],
         sprintf("mode %d: %d of %d runs (representative run %s)%s", mk$mode[i], mk$n_runs[i], sel$n_runs,
                 mk$rep_run[i], if (!is.na(mk$avg_sim_within[i])) sprintf(", mean similarity %.3f", mk$avg_sim_within[i]) else ""))
}))
save_both(bar_plot(x, sprintf("%s: ADMIXTURE modes at K = %d (pong)", opt$species, kcv),
                   paste0("replicate runs grouped by pong; clusters matched to mode 1",
                          if (nrow(mk_all) > MAX_MODES) sprintf("; largest %d of %d modes (%d runs in the others)",
                                                                MAX_MODES, nrow(mk_all), sum(mk_all$n_runs[-seq_len(MAX_MODES)]))
                          else "")),
          "selected_k_modes", 1.6 + 1.6 * nrow(mk))

# CV error
ksum[, stab := factor(ifelse(stable == "yes", "stable", "not stable"), levels = c("stable", "not stable"))]
p <- ggplot(ksum, aes(K, cv_mean)) +
  geom_vline(xintercept = kcv, linetype = "dashed", colour = "grey50") +
  geom_point(data = cv, aes(K, cv_error), inherit.aes = FALSE, colour = "grey60", size = 1.2,
             position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  geom_line() +
  geom_point(aes(shape = stab), size = 3, fill = "white", stroke = 0.9) +
  scale_shape_manual(values = c(stable = 16, `not stable` = 24), name = "pong", drop = FALSE) +
  scale_x_continuous(breaks = ks) +
  labs(x = "K", y = "Cross-validation error",
       title = sprintf("%s: ADMIXTURE cross-validation", opt$species),
       subtitle = sprintf("grey: each of %d runs per K; black: mean\ndashed line: lowest mean CV error (K = %d)",
                          sel$n_runs, kcv))
ggsave(file.path(d, "cv_error.pdf"), p, width = 6.5, height = 4.5)
ggsave(file.path(d, "cv_error.png"), p, width = 6.5, height = 4.5, dpi = 300)
