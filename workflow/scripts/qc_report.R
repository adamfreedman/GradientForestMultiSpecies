#!/usr/bin/env Rscript
# Per-species QC report (self-contained HTML; figures are linked by relative
# path, so keep qc.html inside the species results directory).
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))
suppressPackageStartupMessages(library(ggplot2))

opt <- parse_args()
need(opt, c("species", "spdir", "validation", "snp_filter", "pop_counts", "env_sites",
            "predictors", "heatmap", "tuning", "ntree", "summary", "pvals", "out", "out_tuning_png"))
sp <- opt$species
# path relative to the report's directory (files may live in a parent directory)
rel <- function(path) {
  from <- strsplit(normalizePath(opt$spdir), "/", fixed = TRUE)[[1]]
  to <- strsplit(normalizePath(path), "/", fixed = TRUE)[[1]]
  n <- 0
  while (n < min(length(from), length(to)) && from[n + 1] == to[n + 1]) n <- n + 1
  paste(c(rep("..", length(from) - n), to[-seq_len(n)]), collapse = "/")
}
is_refit <- !is.null(opt$refit) && opt$refit == "1"

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
html_table <- function(dt, digits = 4) {
  dt <- as.data.frame(dt)
  for (j in seq_along(dt)) if (is.double(dt[[j]])) dt[[j]] <- signif(dt[[j]], digits)
  head <- paste0("<tr>", paste0("<th>", esc(names(dt)), "</th>", collapse = ""), "</tr>")
  body <- apply(dt, 1, function(r) paste0("<tr>", paste0("<td>", esc(ifelse(is.na(r), "", trimws(r))), "</td>",
                                                         collapse = ""), "</tr>"))
  paste0("<table>", head, paste(body, collapse = ""), "</table>")
}
img <- function(path, width = "100%") sprintf('<img src="%s" style="max-width:%s">', esc(rel(path)), width)

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
preds <- fread(opt$predictors)

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
    theme_bw(base_size = 10)
  ggsave(opt$out_tuning_png, p, width = 11, height = 3.4, dpi = 150)
  tun_png <- opt$out_tuning_png
} else {
  png(opt$out_tuning_png, width = 400, height = 100); plot.new(); text(0.5, 0.5, "fixed ntree"); dev.off()
}

summ <- fread(opt$summary)
pv <- fread(opt$pvals)
if (pv$n_random[1] < 100) warnings <- c(warnings, sprintf("Only %d randomized fits: p-values are coarse (min %.3g).",
                                                          pv$n_random[1], 1 / (pv$n_random[1] + 1)))

# ---- assemble ----------------------------------------------------------------
plots <- file.path(opt$spdir, "plots")
section <- function(title, ...) paste0("<h2>", esc(title), "</h2>\n", paste(..., sep = "\n"))
pre <- function(x) paste0("<pre>", esc(paste(x, collapse = "\n")), "</pre>")

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
</style></head><body>
<h1>', esc(sp), ': gradient forest QC</h1><p>generated ', format(Sys.time(), "%Y-%m-%d %H:%M"), '</p>',
warn_html,
if (is_refit) paste0('<div class="warn"><b>Refit on significant predictors.</b> These predictors were chosen ',
  'because their importance beat the randomized fits of the full model (', esc(paste(fread(opt$predictors)[kept == TRUE, variable], collapse = ", ")),
  '). The randomization p-values below are conditional on that choice, which used the same data, so they ',
  'are optimistic; the full model (<a href="../qc.html">../qc.html</a>) is the primary test. Use this refit for ',
  'maps, turnover space and per-SNP R&sup2; without the non-significant predictors.</div>'),
section("Model summary", html_table(summ)),
section("Number of trees", ntree_txt, html_table(tun), if (!is.null(tun_png)) img(tun_png)),
section("Observed vs. randomized fits", html_table(pv), img(file.path(plots, "random_vs_real.png")),
        img(file.path(plots, "importance_r2weighted.png"), "49%"), img(file.path(plots, "importance_accuracy.png"), "49%")),
section("Turnover", img(file.path(plots, "map_turnover_rgb.png"), "49%"), img(file.path(plots, "biplot_turnover_pca.png"), "49%")),
section(if (is_refit) "Predictors (selected from the full model)" else "Predictors",
        if (is_refit) "<p>Kept = empirical p &le; alpha for at least one importance metric in the full model. The correlation heat map is from the full predictor set.</p>",
        html_table(preds), img(opt$heatmap, "60%")),
section("SNP filtering", html_table(snpf)),
section("Populations", html_table(merge(sites[, .(popid, sitename, lat, long, snapped_km)], pops, by = "popid", all.x = TRUE)[order(as.numeric(popid))])),
section("Input validation", pre(val)),
"</body></html>")
writeLines(html, opt$out)
if (length(warnings)) message("QC warnings:\n", paste(" -", warnings, collapse = "\n"))
