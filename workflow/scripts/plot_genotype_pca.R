#!/usr/bin/env Rscript
# Genotype PCA plots (PC1 x PC2 and PC1 x PC3): each population has its own
# combination of colour and symbol (legend), with each population's ID at the median position of its individuals (so a few outlying individuals do not move it).
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

opt <- parse_args()
need(opt, c("species", "scores", "variance", "info", "check", "outdir"))
PAIRS <- list(c("PC1", "PC2"), c("PC1", "PC3"))

theme_set(theme_bw(base_size = 11) +  # black panel border, no gridlines
            theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA)))

sc <- fread(opt$scores, colClasses = list(character = c("sample_id", "popid", "sitename")))
va <- fread(opt$variance)
info <- fread(opt$info)
n_snps <- info[item == "snps_used", n]
# with very few individuals there may be fewer than 3 PCs: plot a missing PC as 0
for (pc in setdiff(unique(unlist(PAIRS)), names(sc))) {
  sc[, (pc) := 0]
  va <- rbind(va, data.table(PC = pc, eigenvalue = 0, pct_variance = 0, cumulative_pct = 100))
}

# numeric popids sort as numbers
pop_levels <- unique(sc$popid)
num <- suppressWarnings(as.numeric(pop_levels))
pop_levels <- if (all(!is.na(num))) pop_levels[order(num)] else sort(pop_levels)
sc[, popid := factor(popid, levels = pop_levels)]

# colour x symbol: 8 colours (Okabe-Ito, without yellow, which is hard to see on white)
# and 7 symbols. Population i gets colour i mod 8 and symbol i mod 7, so neighbouring
# IDs differ in both, and the first 56 populations each have a unique combination.
cols <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#999999")
shps <- c(16, 17, 15, 18, 4, 8, 3)
n_pop <- length(pop_levels)
idx <- seq_len(n_pop) - 1
pop_col <- setNames(cols[idx %% length(cols) + 1], pop_levels)
pop_shp <- setNames(shps[idx %% length(shps) + 1], pop_levels)
if (n_pop > length(cols) * length(shps))
  message("WARNING: more than ", length(cols) * length(shps), " populations; colour/symbol combinations repeat")
leg_ncol <- ceiling(n_pop / 22)
leg_w <- leg_ncol * (0.45 + 0.075 * max(nchar(pop_levels), 10))

chk <- fread(opt$check, colClasses = list(character = c("sample_id", "popid")))
flagged <- sc[sample_id %in% chk[flagged == "yes", sample_id]]

axis_lab <- function(pc) sprintf("%s (%.1f%% of variance)", pc, va[PC == pc, pct_variance])
for (pr in PAIRS) {
  x <- pr[1]; y <- pr[2]
  cent <- sc[, .(cx = median(get(x)), cy = median(get(y))), by = popid]
  p <- ggplot(sc, aes(.data[[x]], .data[[y]], colour = popid, shape = popid)) +
    geom_point(size = 2.2, alpha = 0.8) +
    geom_point(data = flagged, shape = 21, size = 5, stroke = 0.8, colour = "black", fill = NA,
               show.legend = FALSE) +
    geom_label_repel(data = cent, aes(cx, cy, label = popid), inherit.aes = FALSE, colour = pop_col[as.character(cent$popid)],
                     size = 3.5, fontface = "bold", fill = alpha("white", 0.8), label.size = 0.25,
                     box.padding = 0.3, min.segment.length = 0, seed = 1, max.overlaps = Inf, show.legend = FALSE) +
    scale_colour_manual(name = "Population", values = pop_col, drop = FALSE) +
    scale_shape_manual(name = "Population", values = pop_shp, drop = FALSE) +
    guides(colour = guide_legend(ncol = leg_ncol, override.aes = list(size = 2.6, alpha = 1)),
           shape = guide_legend(ncol = leg_ncol)) +
    labs(x = axis_lab(x), y = axis_lab(y),
         title = sprintf("%s: genotype PCA", opt$species),
         subtitle = paste0(sprintf("%d individuals, %s SNPs; labels = population ID at the population median",
                                   nrow(sc), format(n_snps, big.mark = ",")),
                           if (nrow(flagged)) sprintf("\nblack rings: %d individual(s) flagged as closer to another population (see qc.html)",
                                                      nrow(flagged)) else "")) +
    theme(legend.key.height = unit(0.9, "lines"))
  out <- file.path(opt$outdir, sprintf("pca_%s_%s", tolower(x), tolower(y)))
  ggsave(paste0(out, ".pdf"), p, width = 7.5 + leg_w, height = 6.5)
  ggsave(paste0(out, ".png"), p, width = 7.5 + leg_w, height = 6.5, dpi = 300)
}
