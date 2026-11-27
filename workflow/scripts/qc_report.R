#!/usr/bin/env Rscript
# Per-species QC report: one self-contained HTML file (figures are embedded as
# base64 PNG data, so qc.html can be copied or emailed on its own).
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))
suppressPackageStartupMessages(library(ggplot2))

opt <- parse_args()
need(opt, c("species", "spdir", "validation", "snp_filter", "pop_counts", "env_sites",
            "predictors", "heatmap", "tuning", "ntree", "summary", "pvals", "snp_pvals", "unmapped",
            "out", "out_tuning_png"))
sp <- opt$species
is_refit <- !is.null(opt$refit) && opt$refit == "1"

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
html_table <- function(dt, digits = 4) {
  dt <- as.data.frame(dt)
  # plain numbers (no scientific notation unless very large or very small), per cell
  fmt <- function(v) vapply(v, function(x) {
    if (is.na(x)) return("")
    x <- signif(x, digits)
    format(x, scientific = x != 0 && (abs(x) >= 1e7 || abs(x) < 1e-5), trim = TRUE, big.mark = "")
  }, character(1))
  for (j in seq_along(dt)) if (is.double(dt[[j]])) dt[[j]] <- fmt(dt[[j]])
  head <- paste0("<tr>", paste0("<th>", esc(names(dt)), "</th>", collapse = ""), "</tr>")
  body <- apply(dt, 1, function(r) paste0("<tr>", paste0("<td>", esc(ifelse(is.na(r), "", trimws(r))), "</td>",
                                                         collapse = ""), "</tr>"))
  paste0("<table>", head, paste(body, collapse = ""), "</table>")
}
# base64 encoding in base R (avoids adding a package to the conda env)
b64 <- function(path) {
  raw <- readBin(path, "raw", file.info(path)$size)
  tab <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", "")[[1]]
  n <- length(raw)
  pad <- (3 - n %% 3) %% 3
  v <- as.integer(c(raw, as.raw(rep(0L, pad))))
  m <- matrix(v, nrow = 3)
  w <- m[1, ] * 65536L + m[2, ] * 256L + m[3, ]
  idx <- rbind(w %/% 262144L, (w %/% 4096L) %% 64L, (w %/% 64L) %% 64L, w %% 64L)
  ch <- tab[as.vector(idx) + 1L]
  if (pad) ch[(length(ch) - pad + 1):length(ch)] <- "="
  paste(ch, collapse = "")
}
img <- function(path, width = "100%") {
  if (!file.exists(path)) return(sprintf("<p class='bad'>missing figure: %s</p>", esc(basename(path))))
  sprintf('<img src="data:image/png;base64,%s" alt="%s" style="max-width:%s">', b64(path), esc(basename(path)), width)
}

warnings <- character(0)

# ---- inputs ------------------------------------------------------------------
val <- readLines(opt$validation)
w_start <- grep("^# warnings", val)
e_start <- grep("^# errors", val)
val_warn <- if (length(w_start) && length(e_start)) val[(w_start + 1):(e_start - 1)] else character(0)
val_warn <- val_warn[nzchar(val_warn)]
warnings <- c(warnings, if (length(val_warn)) paste("Input validation:", val_warn))

pops <- fread(opt$pop_counts, colClasses = list(character = "popid"))
sites <- fread(opt$env_sites, colClasses = list(character = "popid"))
dropped <- setdiff(sites$popid, pops$popid)
if (length(dropped)) warnings <- c(warnings, sprintf("Populations dropped for too few individuals: %s",
                                                     paste(dropped, collapse = ", ")))
snpf <- fread(opt$snp_filter)
unm <- fread(opt$unmapped, colClasses = "character")
snp_p <- fread(opt$snp_pvals)
snp_n <- snp_p$n_random[1]
snp_tab <- data.table(
  measure = c("SNPs in the input (popfreq.tsv)", "SNPs tested (observed R² > 0)",
              sprintf("p ≤ %s", c(0.01, 0.05)), sprintf("BH-adjusted p ≤ %s", c(0.01, 0.05, 0.1))),
  n_snps = c(nrow(snp_p), sum(snp_p$rsq > 0),
             sum(snp_p$p_empirical <= 0.01, na.rm = TRUE), sum(snp_p$p_empirical <= 0.05, na.rm = TRUE),
             sum(snp_p$p_bh <= 0.01, na.rm = TRUE), sum(snp_p$p_bh <= 0.05, na.rm = TRUE),
             sum(snp_p$p_bh <= 0.1, na.rm = TRUE)))
preds <- fread(opt$predictors)
if (!is_refit) preds <- preds[order(status != "kept", rank)]

# ---- ntree selection ---------------------------------------------------------
tun <- fread(opt$tuning)
nsel <- fread(opt$ntree)
if (nzchar(nsel$warning[1]) && !is.na(nsel$warning[1])) warnings <- c(warnings, nsel$warning[1])
tun_png <- NULL
if (identical(nsel$mode[1], "auto") && nrow(tun) > 1) {
  thr <- list(snp_rsq_cor = as.numeric(opt$min_snp_rsq_cor),
              importance_max_change = as.numeric(opt$max_importance_change),
              turnover_dist_cor = as.numeric(opt$min_turnover_cor))
  lab <- c(snp_rsq_cor = "per-SNP R² correlation (≥ threshold)",
           importance_max_change = "max importance change (≤ threshold)",
           turnover_dist_cor = "turnover-distance correlation (≥ threshold)")
  long <- melt(tun[!is.na(snp_rsq_cor), .(ntree, snp_rsq_cor, importance_max_change, turnover_dist_cor)],
               id.vars = "ntree")
  long[, panel := factor(lab[as.character(variable)], levels = lab)]
  th <- data.table(panel = factor(lab[names(thr)], levels = lab), y = unlist(thr))
  p <- ggplot(long, aes(ntree, value)) +
    geom_hline(data = th, aes(yintercept = y), colour = "red3", linetype = 2) +
    geom_line() + geom_point(size = 2) +
    geom_vline(xintercept = nsel$ntree[1], colour = "steelblue4", alpha = 0.5, linewidth = 1) +
    scale_x_log10(breaks = tun$ntree) +
    facet_wrap(~panel, scales = "free_y", nrow = 1) +
    labs(x = "ntree (each point compares this fit with the previous, smaller one)", y = NULL,
         title = sprintf("%s: ntree selection (chosen: %d%s)", sp, nsel$ntree[1],
                         if (isTRUE(as.logical(nsel$converged[1]))) "" else ", NOT converged")) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA))
  ggsave(opt$out_tuning_png, p, width = 11, height = 3.4, dpi = 150)
  tun_png <- opt$out_tuning_png
} else {
  png(opt$out_tuning_png, width = 400, height = 100); plot.new(); text(0.5, 0.5, "fixed ntree"); dev.off()
}

summ <- fread(opt$summary)
# median R² of the SNPs in the model (R² > 0), next to the mean
summ[, median_rsq := median(snp_p$rsq[snp_p$rsq > 0 & !is.na(snp_p$rsq)])]
if ("mean_rsq" %in% names(summ)) setcolorder(summ, "median_rsq", after = "mean_rsq")
pv <- fread(opt$pvals)
if (pv$n_random[1] < 100) warnings <- c(warnings, sprintf("Only %d randomized fits: p-values are coarse (min %.3g).",
                                                          pv$n_random[1], 1 / (pv$n_random[1] + 1)))

# ---- assemble ----------------------------------------------------------------
plots <- file.path(opt$spdir, "plots")
section <- function(title, ...) paste0("<h2>", esc(title), "</h2>\n", paste(..., sep = "\n"))
pre <- function(x) paste0("<pre>", esc(paste(x, collapse = "\n")), "</pre>")

# plain-language explanations shown before each figure and table
legend <- function(...) paste0('<div class="legend">', paste0(..., collapse = ""), "</div>")
tech <- function(...) paste0('<p class="tech"><b>Technical detail:</b> ', paste0(..., collapse = ""), "</p>")
n_rand <- pv$n_random[1]
n_snps_in <- summ$n_snps[1]
pv_show <- copy(pv)
pv_show[, meaning := fcase(
  statistic == "n_pos_rsq", "Number of SNPs whose allele frequencies the environment explains at all (R\u00b2 > 0)",
  statistic == "mean_rsq", "Average strength of fit (R\u00b2) for those SNPs",
  statistic == "sum_rsq", "Total explained variation, summed over those SNPs",
  startsWith(statistic, "imp_"), paste0("Importance of ", sub("^imp_", "", statistic), " (R\u00b2-weighted)"),
  startsWith(statistic, "acc_"), paste0("Importance of ", sub("^acc_", "", statistic), " (accuracy)"),
  default = "")]
setcolorder(pv_show, c("statistic", "meaning"))

intro_html <- legend(
  "<p><b>How to read this report.</b> Gradient forest asks whether differences in allele frequencies among ",
  "populations line up with differences in their environment, and which environmental variables matter most. ",
  "For each SNP, it fits a model predicting the SNP's allele frequency in each population from the environmental ",
  "variables at that population's site. <b>R&sup2;</b> is the fraction of the among-population variation in that ",
  "SNP's allele frequency the model explains (0 = none, 1 = all). SNPs whose R&sup2; is above 0 are said to ",
  "<i>enter the model</i>; the rest are ignored.</p>",
  "<p>Because a flexible model can find patterns in random noise, every result here is compared with ",
  sprintf("<b>%d randomized fits</b>: the same analysis repeated after shuffling which population gets which ", n_rand),
  "environment. Shuffling keeps the environmental data and the genetic data each intact but breaks any real link ",
  "between them. A result that the real data beat in (nearly) all randomized fits is unlikely to be chance. ",
  sprintf("The <b>p-value</b> is the fraction of randomized fits that did at least as well as the real data; with %d fits ", n_rand),
  sprintf("the smallest possible p is %.2g.</p>", 1 / (n_rand + 1)))

summary_legend <- legend(
  "<p>The size of the analysis and how well the model fits overall.</p><ul>",
  "<li><b>n_sites</b>: populations used. <b>n_predictors</b>: environmental variables used.</li>",
  "<li><b>n_snps</b>: SNPs given to gradient forest. <b>n_pos_rsq</b>: SNPs that entered the model (R&sup2; &gt; 0). ",
  "<b>mean_rsq</b>, <b>median_rsq</b> and <b>sum_rsq</b>: average, median and total R&sup2; of those SNPs. ",
  "The mean is pulled down when many weakly associated SNPs enter the model, so read it together with n_pos_rsq; ",
  "the median is less affected by them. Only the mean is tested against the randomized fits (table below).</li>",
  "<li><b>n_snps_few_unique</b>: SNPs with 5 or fewer distinct allele frequencies across populations; the software ",
  "warns about these but still uses them. <b>n_snps_constant</b>: SNPs with the same frequency everywhere, which can ",
  "never enter the model. <b>n_rf_errors</b> / <b>n_rf_warnings</b>: problems reported while fitting (details in the gf_fit log).</li>",
  "<li><b>ntree</b>, <b>maxLevel</b>, <b>corr_threshold</b>: model settings (see Number of trees below and config.yaml).</li></ul>")

ntree_legend <- legend(
  "<p>Each SNP's model is an average over many decision trees; more trees give more stable results but take longer. ",
  "The workflow starts with a small number of trees and keeps doubling it until two successive runs give ",
  "essentially the same answer. The table lists each run; the figure shows how closely each run agreed with the ",
  "previous, smaller one.</p>",
  "<p><b>Figure:</b> each panel is one agreement measure, plotted against the number of trees. The red dashed line is ",
  "the agreement required; the blue vertical line marks the number of trees chosen. Left: correlation of per-SNP ",
  "R&sup2; between runs (higher is better). Middle: largest change in any variable's importance (lower is better). ",
  "Right: correlation of the predicted genetic differences between map locations (higher is better).</p>")

random_table_legend <- legend(
  "<p>This table tests the whole model, and each environmental variable, against the randomized fits. ",
  "Each row is one statistic: <b>observed</b> is its value for the real data, <b>null_mean</b> its average over the ",
  sprintf("%d randomized fits, <b>null_q95</b> the value that only 5%% of randomized fits exceeded, and ", n_rand),
  "<b>p_empirical</b> the fraction of randomized fits that matched or beat the real data. Small p (for example ",
  "&le; 0.05) means the real data show more environmental association than chance; an observed value above ",
  "null_q95 corresponds to p &le; 0.05.</p><ul>",
  "<li><b>n_pos_rsq</b>: how many SNPs show any association with the environment (R&sup2; &gt; 0). Random data ",
  "also produce some, so the question is whether the real data produce more.</li>",
  "<li><b>mean_rsq</b>: how strong the association is, on average, for those SNPs.</li>",
  "<li><b>sum_rsq</b>: the total association across all SNPs (roughly n_pos_rsq &times; mean_rsq); combines how many ",
  "SNPs respond and how strongly.</li>",
  "<li><b>imp_<i>variable</i></b> (R&sup2;-weighted importance): how much of the explained variation is due to this ",
  "variable. Each SNP's R&sup2; is divided among the variables according to how much each one helped predict that SNP, ",
  "then averaged over the SNPs in the model, so the imp_ values of all variables add up to mean_rsq. A variable ",
  "with a significant imp_ value explains more genetic variation than it would by chance.</li>",
  "<li><b>acc_<i>variable</i></b> (accuracy importance): how much worse the predictions get, on average over the ",
  "SNPs in the model, when this variable's values are scrambled. It measures the same thing as imp_ but is not ",
  "scaled by each SNP's R&sup2;, so it can be negative (scrambling did not hurt) and its units are those of squared ",
  "allele-frequency error. Agreement between imp_ and acc_ makes a variable's importance more convincing.</li></ul>",
  tech("each randomized fit permutes whole rows of the population-by-environment table (keeping the correlations ",
       "among variables) and refits with the same variables, settings and number of trees as the real model. ",
       "p = (1 + no. randomized &ge; observed) / (n_random + 1), one-sided. imp_ = gradientForest ",
       "<code>importance(type = \"Weighted\")</code>; acc_ = <code>importance(type = \"Accuracy\")</code>, the mean ",
       "randomForest %IncMSE over SNPs with R&sup2; &gt; 0."))

random_fig_legend <- legend(
  sprintf("<p><b>Figure: real model vs. randomized fits.</b> Each panel is one whole-model statistic. The grey histogram shows its values in the %d randomized fits; ", n_rand),
  "the red line is the value for the real data, with its p-value. A red line far to the right of the histogram ",
  "means the real data show much more environmental association than shuffled data. Panels: number of SNPs with ",
  "R&sup2; &gt; 0, their mean R&sup2;, and the total R&sup2;.</p>")

importance_fig_legend <- legend(
  "<p><b>Figures: importance of each environmental variable.</b> Left: R&sup2;-weighted importance; right: ",
  "accuracy importance (both explained in the table above). Each bar is the variable's importance in the real data, ",
  "ordered from most to least important. The thin red line under each bar shows the range of importance in the ",
  "randomized fits (5th to 95th percentile; dot = average). A bar that extends past the right end of its red ",
  "line is more important than expected by chance (dark blue bars, p &le; 0.05); grey bars cannot be distinguished ",
  "from chance.</p>")

turnover_legend <- legend(
  "<p><b>What turnover means.</b> The fitted model predicts, for every location on the map, the genetic makeup ",
  "expected there given its environment. <i>Turnover</i> is how much that predicted genetic makeup changes between ",
  "places: locations with similar predicted genetic makeup are close in <i>turnover space</i>, and locations with ",
  "different makeup are far apart. These are predictions from environment alone, including at unsampled locations.</p>",
  "<p><b>Left, map:</b> each grid cell is coloured by its predicted genetic makeup. Cells of similar colour are ",
  "expected to be genetically similar; sharp colour changes mark where the environment predicts rapid genetic change. ",
  "Colours have no meaning on their own beyond similarity. White circles are the sampled populations; lines are the ",
  "study-area outline and country borders.</p>",
  "<p><b>Right, turnover space:</b> the same grid cells, now placed by their predicted genetic makeup instead of ",
  "their geography (first two principal components, PC1 and PC2; percentages are the share of variation each ",
  "shows). Points are coloured as in the map. Arrows show the environmental variables: an arrow points toward ",
  "cells with high values of that variable, and longer arrows mean the variable drives more of the genetic differences ",
  "shown. Labelled circles are the sampled populations.</p>",
  tech("each predictor is transformed with <code>predict(gf)</code> onto its cumulative-importance scale; a PCA ",
       "(unscaled) on the transformed predictors listed in the map subtitle gives PC1-3, which are mapped to colour as ",
       "red = PC1 + PC2, green = &minus;PC2, blue = PC3 + PC2 &minus; PC1, each rescaled to 0-255 across grid cells ",
       "(Fitzpatrick &amp; Keller 2015). Values beyond the sampled environmental range are extrapolated."))

predictors_legend <- legend(
  "<p>Environmental variables that are strongly correlated with each other carry largely the same information and ",
  "make it hard to tell which one matters, so the workflow keeps only one variable from each strongly correlated ",
  "set: no two kept variables have a correlation (|r|) above the threshold across the sampled sites. ",
  "To decide which one to keep, the model is first fitted with <i>all</i> the variables (a \"screening\" fit), and ",
  "the variables are ranked by how much of the genetic variation each one explains. Going down that ranking, ",
  "each variable is kept unless it is strongly correlated with one already kept. So of two correlated variables, ",
  "the one that fits the genetic data better is kept.</p>",
  "<p>The randomized fits use the same kept variables, so every test in this report asks whether shuffled data ",
  "could produce the fit of <i>this</i> model. Because the variables were chosen using the real data, the ",
  "whole-model and variable-importance p-values are somewhat optimistic; the per-SNP tests are conditional on the ",
  "chosen model by design.</p>",
  "<ul><li><b>status</b>: kept or dropped. <b>rank</b>: position in the ranking (1 = explains the most).</li>",
  "<li><b>screen_importance</b>: the variable's importance in the screening fit with all variables (its share of ",
  "the explained genetic variation, as imp_ above).</li>",
  "<li><b>mean_abs_r</b>: its average correlation with the other variables (0 to 1).</li>",
  "<li><b>reason</b>: for a dropped variable, the kept variable it is too strongly correlated with.</li></ul>",
  "<p><b>Figure:</b> correlation between every pair of variables across the sampled sites (red = positive, blue = ",
  "negative, darker = stronger, value printed in each cell); kept variables are in bold.</p>",
  tech("ranking and threshold are set by <code>predictors.rank_by</code> and <code>predictors.corr_threshold</code> ",
       "in config.yaml; with <code>rank_by: correlation</code> the ranking ignores the genetic data (least correlated ",
       "with the others first)."))

# optional parts of the predictor legend: share column (newer report.tsv) and importance floor
if ("screen_share" %in% names(preds))
  predictors_legend <- paste0(predictors_legend, legend(
    "<p><b>screen_share</b>: the variable's share of the summed screening importance of all variables ",
    "(0.05 = 5% of the explained genetic variation).</p>"))
if (!is.null(opt$importance_floor) && !opt$importance_floor %in% c("", "none")) {
  fl <- as.numeric(opt$importance_floor)
  predictors_legend <- paste0(predictors_legend, legend(sprintf(paste0(
    "<p><b>Importance floor: %g%%</b> (<code>--importance-floor %g</code>). Variables whose screen_share is below ",
    "%g%% were dropped even if not correlated with a kept variable, so that variables explaining almost none of ",
    "the genetic variation do not add noise to the model and the randomization tests (reason column).</p>"),
    100 * fl, fl, 100 * fl)))
}

no_selection_legend <- legend(
  "<p><b>Variable selection was turned off</b> (<code>--use-all-predictors</code>, or ",
  "<code>predictors.corr_threshold: null</code>): every candidate variable that varies across the sampled sites is ",
  "in the model, including strongly correlated ones (see the figure). Importance is then shared among correlated ",
  "variables, so a low importance for one of them does not mean its signal is unimportant.</p>",
  "<p><b>Figure:</b> correlation between every pair of variables across the sampled sites (red = positive, blue = ",
  "negative, darker = stronger, value printed in each cell).</p>")

stab_html <- ""
if (!is_refit && !is.null(opt$stability_vars) && nzchar(opt$stability_vars)) {
  sv <- fread(opt$stability_vars)
  ss <- fread(opt$stability_sites, colClasses = list(character = "left_out"))
  n_fit <- nrow(ss)
  n_same <- sum(ss$same_as_full == "yes")
  kept <- sv[status == "kept"]
  unstable <- kept[kept_in < 0.75 * n_fit, variable]
  if (length(unstable)) warnings <- c(warnings, sprintf(
    "Predictor choice is sensitive to the sampled populations: %s kept in fewer than 75%% of leave-one-out selections (see Predictors).",
    paste(unstable, collapse = ", ")))
  stab_legend <- legend(
    "<p><b>How much does the choice of variables depend on which populations were sampled?</b> ",
    sprintf("The selection above was repeated %d times, each time leaving out one population ", n_fit),
    "(screening fit and correlations both recomputed without it). ",
    sprintf("The full-data selection was reproduced exactly in <b>%d of %d</b> of these.</p>", n_same, n_fit),
    "<p>A kept variable that stays kept in nearly every case is a robust choice. A kept variable that is often ",
    "replaced is best read as standing for its correlated group (listed in the last column) rather than as the ",
    "single driver: which member of the group is chosen depends on particular populations. A variable can also ",
    "drop out with no substitute, when without one population it becomes too correlated with another kept variable. ",
    "This table describes sensitivity only; the model and tests above use the full-data selection.</p><ul>",
    "<li><b>status</b>, <b>full_rank</b>: the variable's status and rank in the full-data selection above.</li>",
    sprintf("<li><b>kept_in</b> / <b>kept_pct</b>: in how many of the %d leave-one-out selections it was kept.</li>", n_fit),
    "<li><b>median_rank</b>, <b>min_rank</b>, <b>max_rank</b>: its rank across those selections.</li>",
    "<li><b>when_not_kept_replaced_by</b>: when it was not kept, the kept variable most correlated with it ",
    "(number of selections; average |r|), i.e. what stood in for it.</li></ul>")
  sites_legend <- legend(
    "<p><b>Per left-out population.</b> <b>n_pos_rsq</b>: SNPs with R&sup2; &gt; 0 in the screening fit without ",
    "that population (a much higher or lower value than for the others marks a population with strong influence ",
    "on the fit). <b>n_kept</b>: variables selected. <b>added</b> / <b>lost</b>: variables selected that are not in ",
    "the full-data selection, and full-data variables not selected. Populations that change the selection most ",
    "are the ones it depends on.</p>",
    tech("set <code>predictors.stability: false</code> in config.yaml to skip this diagnostic (one screening fit ",
         "per population)."))
  stab_html <- paste0("<h3>Sensitivity to the sampled populations (leave one population out)</h3>",
                      stab_legend, html_table(sv), sites_legend, html_table(ss))
}

snp_filter_legend <- legend(
  "<p>SNPs removed before the analysis. Only SNPs that cannot be analysed are removed: those with no genotype data ",
  "in some population (so no allele frequency there) and those with the same allele frequency in every population ",
  "(nothing to explain). There is no filter on allele frequency or missing data; SNPs with too little information ",
  "simply do not enter the model.</p>")

snp_tests_legend <- legend(
  "<p>Which individual SNPs are associated with the environment more strongly than chance. Only SNPs that entered ",
  "the model (R&sup2; &gt; 0) are tested. A SNP's p-value is the fraction of randomized fits in which that SNP got an ",
  "R&sup2; at least as high as in the real data. Because many SNPs are tested at once, some small p-values are ",
  "expected by chance; the <b>BH-adjusted p</b> (Benjamini-Hochberg) corrects for this, so that among SNPs with BH ",
  "&le; 0.05, about 5% are expected to be false positives. The full per-SNP table, with positions, is ",
  "<code>gf/snp_pvalues.tsv</code>.</p>",
  tech(sprintf("p = (1 + no. of the %d randomized fits with R&sup2; &ge; observed) / (%d + 1); smallest possible p = %.2g. ",
               snp_n, snp_n, 1 / (snp_n + 1)),
       "BH is applied over the tested SNPs only; SNPs with R&sup2; &le; 0 have no p-value."))

pops_legend <- legend("<p>Populations in the analysis: location, number of individuals with genotypes (n_ind), ",
                      "and snapped_km, the distance a site was moved to the nearest grid cell with environmental data ",
                      "(0 = not moved). Populations with too few individuals are left out (see Warnings).</p>")

# genotype PCA (individuals; qc.html only)
pca_html <- NULL
if (!is.null(opt$pca_png)) {
  pva <- fread(opt$pca_variance)
  psc <- fread(opt$pca_scores, colClasses = list(character = c("sample_id", "popid", "sitename")))
  pinf <- fread(opt$pca_info)
  hi_miss <- psc[missing_frac > 0.5]
  if (nrow(hi_miss)) warnings <- c(warnings, sprintf(
    "Genotype PCA: %d individual(s) with more than 50%% missing genotypes (pulled towards the centre of the PCA): %s",
    nrow(hi_miss), paste(head(hi_miss$sample_id, 10), collapse = ", ")))
  not_gf <- unique(psc[in_gf_model == "no", popid])
  pca_legend <- legend(
    "<p>Principal components of the individual genotypes: a check on population structure and on the samples, ",
    "independent of the environment and of gradient forest. Left: PC1 against PC2; right: PC1 against PC3. ",
    "Each point is an individual; its position comes from its genotypes, and its colour and symbol from its ",
    "population in the popmap (each population has its own combination; see the legend). The population ID is ",
    "placed at the median of its individuals (each PC separately). Individuals of the same population usually ",
    "cluster together. Black rings mark individuals flagged as sitting with another population (table below); ",
    "an individual near the centre with high missing_frac may simply lack data.</p>",
    sprintf("<p>%d individuals and %s SNPs (biallelic SNPs variable among these individuals; no MAF or missingness ",
            pinf[item == "n_individuals", n], format(pinf[item == "snps_used", n], big.mark = ",")),
    "filter). The table gives the share of the genotype variance captured by the first PCs. Per-individual ",
    "coordinates on every PC are in <code>pca/sample_scores.tsv</code>.",
    if (length(not_gf)) sprintf(" Population(s) %s are included here but not in the gradient forest model (too few individuals).",
                                esc(paste(not_gf, collapse = ", "))),
    "</p>",
    tech("SNPs are centred and scaled by sqrt(2p(1-p)) (Patterson et al. 2006); missing genotypes are set to the mean."))
  pchk <- fread(opt$pca_check, colClasses = list(character = c("sample_id", "popid", "nearest_other_pop", "note")))
  cs <- as.numeric(strsplit(opt$pca_check_settings, ",")[[1]])
  flag <- pchk[flagged == "yes"][order(-ratio)]
  if (nrow(flag)) warnings <- c(warnings, sprintf(
    "Genotype PCA: %d individual(s) sit with another population rather than their own (possible popmap or sample errors): %s",
    nrow(flag), paste(head(flag$sample_id, 10), collapse = ", ")))
  pchk[is.na(note), note := ""]
  n_unchecked <- sum(pchk$note %in% "only individual in its population: not checked")
  check_html <- paste0(
    legend(sprintf(paste0(
      "<p><b>Individuals that sit with another population.</b> On the first %d PCs, an individual is flagged ",
      "when it is at least %g times closer to another population's median than to its own (<i>ratio</i>), and ",
      "more than %g times farther from its own population's median than that population's members usually are ",
      "(<i>spread_ratio</i>). Both conditions are needed because similar populations overlap: many individuals ",
      "sit a little closer to a neighbour without being out of place. Medians are computed without the ",
      "individual being checked; for a population of fewer than 3, the typical spread of all populations is used. A flagged individual may be mislabelled in the popmap, a sample mix-up or ",
      "contaminated, a migrant or admixed individual, or from a genetically mixed population. ",
      "All individuals, flagged or not, are in <code>pca/population_check.tsv</code>; the thresholds are ",
      "<code>pca.population_check</code> in config.yaml.%s</p>"), cs[1], cs[2], cs[3],
      if (n_unchecked) sprintf(" %d individual(s) alone in their population could not be checked.", n_unchecked) else "")),
    if (nrow(flag)) html_table(merge(flag[, .(sample_id, popid, nearest_other_pop, ratio, spread_ratio, note)],
                                     psc[, .(sample_id, missing_frac)], by = "sample_id", sort = FALSE), digits = 3)
    else "<p>No individuals flagged.</p>")
  pca_html <- section("Genotype PCA", pca_legend, paste(vapply(strsplit(opt$pca_png, ",")[[1]], img, "", width = "49%"), collapse = " "),
                      html_table(head(pva[, .(PC, pct_variance, cumulative_pct)], 10), digits = 3),
                      check_html,
                      if (nrow(hi_miss)) paste0("<p>Individuals with more than 50% missing genotypes:</p>",
                                                html_table(hi_miss[, .(sample_id, popid, missing_frac)], digits = 3)))
}

warn_html <- if (length(warnings)) {
  paste0('<div class="warn"><b>Warnings</b><ul>', paste0("<li>", esc(warnings), "</li>", collapse = ""), "</ul></div>")
} else '<div class="ok">No warnings.</div>'

ntree_txt <- sprintf("<p>Selected <b>ntree = %d</b> (%s%s). The same ntree is used for all %d randomized fits.</p>",
                     nsel$ntree[1], nsel$mode[1],
                     if (identical(nsel$mode[1], "auto"))
                       if (isTRUE(as.logical(nsel$converged[1]))) ", converged" else
                         sprintf(", <span class='bad'>stopped at maximum %d without converging</span>", nsel$tune_max[1])
                     else "",
                     pv$n_random[1])

html <- paste0('<!doctype html><html><head><meta charset="utf-8"><title>', esc(sp), ' QC</title>
<style>
body{font-family:-apple-system,Helvetica,Arial,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;color:#222}
h1{margin-bottom:0} h2{border-bottom:1px solid #ccc;padding-bottom:4px;margin-top:32px}
table{border-collapse:collapse;font-size:13px;margin:8px 0} th,td{border:1px solid #ddd;padding:3px 8px;text-align:left}
th{background:#f3f3f3} pre{background:#f7f7f7;padding:8px;font-size:12px;overflow-x:auto}
.warn{background:#fff4e5;border-left:4px solid #e69500;padding:8px 12px} .ok{background:#eef8ee;border-left:4px solid #3a3;padding:8px 12px}
.bad{color:#b00;font-weight:bold} img{border:1px solid #eee;margin:6px 0}
.legend{background:#f4f7fb;border-left:4px solid #6b8fb8;padding:4px 12px;margin:8px 0;font-size:14px;line-height:1.45}
.legend ul{margin:4px 0 8px 0;padding-left:20px} .legend li{margin:3px 0}
.tech{font-size:12px;color:#555}
</style></head><body>
<h1>', esc(sp), ': gradient forest QC</h1><p>generated ', format(Sys.time(), "%Y-%m-%d %H:%M"), '</p>',
warn_html,
if (is_refit) paste0('<div class="warn"><b>Refit on significant predictors.</b> These predictors were chosen ',
  'because their importance beat the randomized fits of the full model (', esc(paste(fread(opt$predictors)[kept == TRUE, variable], collapse = ", ")),
  '). The randomization p-values below are conditional on that choice, which used the same data, so they ',
  'are optimistic; the full model (<a href="../qc.html">../qc.html</a>) is the primary test. Use this refit for ',
  'maps, turnover space and per-SNP R&sup2; without the non-significant predictors.</div>'),
intro_html,
section("Model summary", summary_legend, html_table(summ)),
section("Observed vs. randomized fits", random_table_legend, html_table(pv_show),
        random_fig_legend, img(file.path(plots, "random_vs_real.png")),
        importance_fig_legend,
        img(file.path(plots, "importance_r2weighted.png"), "49%"), img(file.path(plots, "importance_accuracy.png"), "49%")),
section("Turnover", turnover_legend,
        img(file.path(plots, "map_turnover_rgb.png"), "49%"), img(file.path(plots, "biplot_turnover_pca.png"), "49%")),
section("Per-SNP tests", snp_tests_legend, html_table(snp_tab)),
section(if (is_refit) "Predictors (selected from the full model)" else "Predictors",
        if (is_refit) legend("<p>In this refit, the variables kept are those whose importance beat the randomized fits ",
                             "of the full model (p &le; alpha for at least one importance measure). The correlation figure ",
                             "is from the full set of variables.</p>")
        else if (identical(opt$selection, "none")) no_selection_legend else predictors_legend,
        html_table(preds), img(opt$heatmap, "60%"), stab_html),
section("SNP filtering", snp_filter_legend, html_table(snpf)),
section("Number of trees", ntree_legend, ntree_txt, html_table(tun), if (!is.null(tun_png)) img(tun_png)),
pca_html,
section("Populations", pops_legend, html_table(merge(sites[, .(popid, sitename, lat, long, snapped_km)], pops, by = "popid", all.x = TRUE)[order(suppressWarnings(as.numeric(popid)), popid)])),
section("Individuals not mapped to a population",
        if (nrow(unm)) paste0("<p>Excluded from the analysis. <i>missing_in_vcf</i> = in the popmap but not the VCF; ",
                              "<i>missing_in_popmap</i> = in the VCF but not the popmap; <i>popid_in_latlongs</i> = ",
                              "whether the individual's population has coordinates.</p>", html_table(unm))
        else "<p>All individuals in the VCF and the popmap map to a population with coordinates.</p>"),
section("Input validation", pre(val)),
"</body></html>")
writeLines(html, opt$out)
if (length(warnings)) message("QC warnings:\n", paste(" -", warnings, collapse = "\n"))
