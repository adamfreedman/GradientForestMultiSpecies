#!/usr/bin/env Rscript
# Single-species gradient forest figures.
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))
suppressPackageStartupMessages({
  library(gradientForest)
  library(terra)
  library(ggplot2)
  library(ggrepel)
  library(grid)
})

opt <- parse_args()
need(opt, c("species", "model", "importance", "summary", "null", "pvals", "grid", "pca",
            "site_pcs", "study_area", "boundaries", "biplot_max_points", "outdir"))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
sp <- opt$species

save_plot <- function(p, name, w, h) {
  ggsave(file.path(opt$outdir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(opt$outdir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
theme_set(theme_bw(base_size = 11) +  # black panel border, no gridlines
            theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA)))

gf <- readRDS(opt$model)
summ <- fread(opt$summary)
imp <- fread(opt$importance)
null <- fread(opt$null)
pv <- fread(opt$pvals)
grid <- as.data.table(readRDS(opt$grid))
pinfo <- readRDS(opt$pca)
sites <- fread(opt$site_pcs, colClasses = list(character = "popid"))
grid[, col := rgb(r, g, b, maxColorValue = 255)]
sites[, col := rgb(r, g, b, maxColorValue = 255)]

# study-area outline as a path data frame (avoids an sf dependency)
area <- as.data.frame(geom(vect(opt$study_area)))
area$grp <- paste(area$geom, area$part, area$hole)

# country boundaries (outlines of the clipped polygons), drawn over the turnover colours
bnd <- tryCatch(vect(opt$boundaries), error = function(e) vect())  # empty layer when boundaries: none
borders <- if (nrow(bnd) > 0) {
  b <- as.data.frame(geom(as.lines(bnd)))
  b$grp <- paste(b$geom, b$part)
  b
} else NULL

# ---- 1. turnover map -------------------------------------------------------
p_map <- ggplot() +
  geom_raster(data = grid, aes(x, y, fill = col)) +
  scale_fill_identity() +
  geom_path(data = area, aes(x, y, group = grp), colour = "grey30", linewidth = 0.3) +
  (if (!is.null(borders)) geom_path(data = borders, aes(x, y, group = grp), colour = "grey10",
                                    linewidth = 0.35, alpha = 0.8)) +
  geom_point(data = sites, aes(long, lat), shape = 21, fill = "white", colour = "black", size = 2.2) +
  geom_text_repel(data = sites, aes(long, lat, label = sitename), size = 3, min.segment.length = 0,
                  bg.colour = "white", bg.r = 0.15, seed = 1) +
  coord_quickmap(xlim = range(grid$x), ylim = range(grid$y), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude",
       title = sprintf("%s: predicted genomic turnover", sp),
       subtitle = paste(strwrap(sprintf("%s grid cells; colours = PCA of transformed %s",
                                        format(nrow(grid), big.mark = ","),
                                        paste(pinfo$pca_vars, collapse = ", ")), 80), collapse = "\n"))
save_plot(p_map, "map_turnover_rgb", 8, 8)

# ---- 2. PCA biplot of turnover space --------------------------------------
set.seed(1)
nmax <- as.integer(opt$biplot_max_points)
gsub_ <- if (nrow(grid) > nmax) grid[sample(.N, nmax)] else grid
load <- as.data.table(pinfo$pca$rotation[, 1:2, drop = FALSE], keep.rownames = "variable")
scal <- 0.9 * max(abs(c(gsub_$PC1, gsub_$PC2))) / max(sqrt(load$PC1^2 + load$PC2^2))
load[, `:=`(xe = PC1 * scal, ye = PC2 * scal)]
ve <- 100 * pinfo$var_explained
p_bi <- ggplot() +
  geom_point(data = gsub_, aes(PC1, PC2), colour = gsub_$col, size = 0.4, alpha = 0.6, shape = 16) +
  geom_segment(data = load, aes(0, 0, xend = xe, yend = ye),
               arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.5) +
  geom_text_repel(data = load, aes(xe * 1.08, ye * 1.08, label = variable), size = 3.2,
                  fontface = "bold", min.segment.length = Inf, seed = 1) +
  geom_point(data = sites, aes(PC1, PC2), shape = 21, fill = sites$col, colour = "black", size = 2.6) +
  geom_text_repel(data = sites, aes(PC1, PC2, label = sitename), size = 2.8, colour = "grey20",
                  bg.colour = "white", bg.r = 0.15, seed = 2) +
  coord_equal() +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = sprintf("%s: turnover space", sp),
       subtitle = "grid cells coloured as in the map; arrows = predictor loadings; circles = sampled populations")
save_plot(p_bi, "biplot_turnover_pca", 8, 7)

# ---- 3. variable importance vs. randomizations -----------------------------
# Bar = observed importance. Bars are coloured by the empirical p-value of the
# observed value against the randomized (environment-permuted) fits; the red
# tick marks the 95th percentile of the randomized values. The smallest
# attainable p is 1 / (n_random + 1).
# colour key drawn inside the panel (bottom right, where the smallest bars are)
importance_key <- function(n) {
  row_y <- function(i) unit(1, "npc") - unit(4 + 5 * i, "mm")
  fs <- gpar(fontsize = 7.5)
  # annotation_custom() overwrites the outer grob's viewport, so the key's own
  # viewport sits on an inner gTree
  key <- gTree(
    vp = viewport(x = unit(1, "npc") - unit(3, "mm"), y = unit(3, "mm"),
                  width = unit(92, "mm"), height = unit(37, "mm"), just = c("right", "bottom")),
    children = gList(
      rectGrob(gp = gpar(fill = "white", col = "grey50", alpha = 0.92)),
      textGrob(sprintf("Bar = observed importance vs. %d randomized fits", n),
               x = unit(3, "mm"), y = row_y(0), just = "left", gp = gpar(fontsize = 7.5, fontface = "bold")),
      rectGrob(x = unit(3, "mm"), y = row_y(1), width = unit(5, "mm"), height = unit(3, "mm"),
               just = c("left", "centre"), gp = gpar(fill = "steelblue4", col = NA)),
      textGrob("higher than randomized fits (empirical p <= 0.05)", x = unit(10, "mm"), y = row_y(1),
               just = "left", gp = fs),
      rectGrob(x = unit(3, "mm"), y = row_y(2), width = unit(5, "mm"), height = unit(3, "mm"),
               just = c("left", "centre"), gp = gpar(fill = "grey65", col = NA)),
      textGrob("not distinguishable from randomized fits (p > 0.05)", x = unit(10, "mm"), y = row_y(2),
               just = "left", gp = fs),
      segmentsGrob(x0 = unit(2.5, "mm"), x1 = unit(8.5, "mm"), y0 = row_y(3), y1 = row_y(3),
                   gp = gpar(col = "red3", lwd = 1.5)),
      pointsGrob(x = unit(5.5, "mm"), y = row_y(3), pch = 16, size = unit(1.8, "mm"),
                 gp = gpar(col = "red3")),
      textGrob("randomized fits: dot = mean, line = 5th to 95th percentile", x = unit(10, "mm"),
               y = row_y(3), just = "left", gp = fs),
      textGrob("bar past the right end of its red line: p <= 0.05 (one-sided)", x = unit(10, "mm"),
               y = row_y(3.8), just = "left", gp = fs),
      textGrob(sprintf("p = (1 + no. randomized >= observed) / (%d + 1)", n),
               x = unit(3, "mm"), y = row_y(5.0), just = "left", gp = gpar(fontsize = 6.5, col = "grey30")),
      textGrob(sprintf("smallest possible p = %.2g", 1 / (n + 1)),
               x = unit(3, "mm"), y = row_y(5.7), just = "left", gp = gpar(fontsize = 6.5, col = "grey30"))
    )
  )
  gTree(children = gList(key))
}

plot_importance <- function(col, prefix, xlab, title, file) {
  d <- imp[, .(variable, value = get(col))]
  st <- pv[startsWith(statistic, prefix)][, variable := sub(prefix, "", statistic, fixed = TRUE)]
  d <- merge(d, st[, .(variable, null_q95, p_empirical)], by = "variable", all.x = TRUE)
  # null mean and 5th percentile from the randomized fits (95th percentile as used for the p-values)
  nc <- paste0(prefix, d$variable)
  d[, null_mean := vapply(nc, function(v) mean(null[[v]], na.rm = TRUE), numeric(1))]
  d[, null_q05 := vapply(nc, function(v) unname(quantile(null[[v]], 0.05, na.rm = TRUE)), numeric(1))]
  d[, variable := factor(variable, levels = variable[order(value)])]
  d[, sig := factor(ifelse(!is.na(p_empirical) & p_empirical <= 0.05,
                           "p <= 0.05", "p > 0.05"), levels = c("p <= 0.05", "p > 0.05"))]
  p <- ggplot(d, aes(value, variable)) +
    geom_col(aes(fill = sig), width = 0.6) +
    # randomized-fit summary drawn just under each bar
    geom_segment(aes(x = null_q05, xend = null_q95, y = as.numeric(variable) - 0.42,
                     yend = as.numeric(variable) - 0.42), colour = "red3", linewidth = 0.6) +
    geom_point(aes(x = null_mean, y = as.numeric(variable) - 0.42), colour = "red3", size = 1.6) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = c("p <= 0.05" = "steelblue4", "p > 0.05" = "grey65"), drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    scale_y_discrete(expand = expansion(add = c(5.1, 0.6))) +  # empty space below the bars for the key
    annotation_custom(importance_key(nrow(null))) +
    labs(x = xlab, y = NULL, title = sprintf("%s: %s", sp, title)) +
    theme(legend.position = "none")
  # leave room under the bars for the key
  save_plot(p, file, 7.5, max(4.5, 0.32 * nrow(d) + 2.6))
}
plot_importance("importance", "imp_", expression(R^2 * "-weighted importance"),
                "R\u00b2-weighted importance", "importance_r2weighted")
plot_importance("accuracy", "acc_", "Accuracy importance (mean over SNPs)",
                "accuracy importance", "importance_accuracy")

# ---- 4. observed vs. randomized model statistics ----------------------------
lab <- c(n_pos_rsq = "Number of SNPs with R^2 > 0", mean_rsq = "Mean R^2 of those SNPs",
         sum_rsq = "Total R^2 (sum over SNPs)")
nl <- melt(null[, names(lab), with = FALSE], measure.vars = names(lab), variable.name = "statistic")
ob <- pv[statistic %in% names(lab)]
nl[, statistic := factor(lab[as.character(statistic)], levels = lab)]
ob[, statistic := factor(lab[statistic], levels = lab)]
ob[, txt := sprintf("observed = %s\np = %.3g", signif(observed, 4), p_empirical)]
p_rand <- ggplot(nl, aes(value)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "grey40") +
  geom_vline(data = ob, aes(xintercept = observed), colour = "red3", linewidth = 0.8) +
  geom_text(data = ob, aes(x = observed, y = Inf, label = txt), hjust = 1.05, vjust = 1.3,
            size = 3, colour = "red3") +
  facet_wrap(~statistic, scales = "free", nrow = 1) +
  expand_limits(x = 0) +
  labs(x = NULL, y = "Randomized fits",
       title = sprintf("%s: observed model vs. %d environment-permuted fits", sp, nrow(null)))
save_plot(p_rand, "random_vs_real", 12, 4)

# ---- 5. standard gradientForest diagnostics ---------------------------------
pdf(file.path(opt$outdir, "gf_diagnostics.pdf"), width = 9, height = 8)
for (pt in c("Overall.Importance", "Split.Density", "Cumulative.Importance", "Performance")) {
  tryCatch({
    if (pt == "Cumulative.Importance") {
      plot(gf, plot.type = pt, show.species = FALSE, common.scale = TRUE,
           imp.vars = names(importance(gf)))
    } else if (pt == "Split.Density") {
      plot(gf, plot.type = pt, imp.vars = names(importance(gf)))
    } else {
      plot(gf, plot.type = pt)
    }
  }, error = function(e) {
    plot.new(); title(main = pt); text(0.5, 0.5, paste("plot failed:", conditionMessage(e)))
  })
}
dev.off()
message("done")
