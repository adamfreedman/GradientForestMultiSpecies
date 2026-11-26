#!/usr/bin/env Rscript
# Genotype PCA plot: individuals on PC1 x PC2, coloured by population, with each
# population's ID at the mean position of its individuals.
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

opt <- parse_args()
need(opt, c("species", "scores", "variance", "info", "out_prefix"))

theme_set(theme_bw(base_size = 11) +  # black panel border, no gridlines
            theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA)))

sc <- fread(opt$scores, colClasses = list(character = c("sample_id", "popid", "sitename")))
va <- fread(opt$variance)
info <- fread(opt$info)
n_snps <- info[item == "snps_used", n]
if (!"PC2" %in% names(sc)) {  # only one PC (e.g. 2 individuals' worth of structure): plot PC1 against 0
  sc[, PC2 := 0]
  va <- rbind(va, data.table(PC = "PC2", eigenvalue = 0, pct_variance = 0, cumulative_pct = 100))
}

# numeric popids sort as numbers
pop_levels <- unique(sc$popid)
num <- suppressWarnings(as.numeric(pop_levels))
pop_levels <- if (all(!is.na(num))) pop_levels[order(num)] else sort(pop_levels)
sc[, popid := factor(popid, levels = pop_levels)]
cent <- sc[, .(PC1 = mean(PC1), PC2 = mean(PC2), n = .N), by = popid]

axis_lab <- function(pc) sprintf("%s (%.1f%% of variance)", pc, va[PC == pc, pct_variance])
p <- ggplot(sc, aes(PC1, PC2, colour = popid)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_label_repel(data = cent, aes(label = popid), size = 3.5, fontface = "bold", fill = alpha("white", 0.8),
                   label.size = 0.25, box.padding = 0.3, min.segment.length = 0, seed = 1,
                   max.overlaps = Inf, show.legend = FALSE) +
  labs(x = axis_lab("PC1"), y = axis_lab("PC2"),
       title = sprintf("%s: genotype PCA", opt$species),
       subtitle = sprintf("%d individuals, %s SNPs; labels = population ID at the population mean",
                          nrow(sc), format(n_snps, big.mark = ","))) +
  theme(legend.position = "none")

ggsave(paste0(opt$out_prefix, ".pdf"), p, width = 7.5, height = 6.5)
ggsave(paste0(opt$out_prefix, ".png"), p, width = 7.5, height = 6.5, dpi = 300)
