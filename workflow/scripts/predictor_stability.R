#!/usr/bin/env Rscript
# Gather the leave-one-population-out predictor selections into two tables:
#   --out_variables: per candidate variable, how often it is kept when one
#     population is left out, its rank range, and what stood in for it;
#   --out_sites: per left-out population, which variables join or leave the
#     selection compared with the full data.
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))

opt <- parse_args()
need(opt, c("selected", "report", "out_variables", "out_sites"))
files <- opt$positional
loo <- rbindlist(lapply(files, fread, colClasses = list(character = c("left_out", "sitename", "variable",
                                                                     "nearest_kept"))), fill = TRUE)
full <- readLines(opt$selected)
full <- full[nzchar(full)]
rep <- fread(opt$report)
n_fit <- uniqueN(loo$left_out)
if (!n_fit) stop("no leave-one-out fits found")

# substitutes: in the fits where a variable is not kept, which kept variable is
# most correlated with it (counted over fits)
subst <- loo[kept == FALSE & !is.na(nearest_kept),
             .(n = .N, r = round(mean(r_nearest), 2)), by = .(variable, nearest_kept)][order(variable, -n)]
subst_txt <- subst[, .(instead = paste(sprintf("%s (%d; |r| %.2f)", nearest_kept, n, r), collapse = "; ")),
                   by = variable]

vars <- loo[, .(kept_in = sum(kept), n_fits = .N,
                median_rank = as.numeric(median(rank, na.rm = TRUE)),
                min_rank = suppressWarnings(min(rank, na.rm = TRUE)),
                max_rank = suppressWarnings(max(rank, na.rm = TRUE))), by = variable]
vars[!is.finite(min_rank), c("min_rank", "max_rank", "median_rank") := NA]
vars <- merge(rep[, .(variable, status, full_rank = rank)], vars, by = "variable", all = TRUE)
vars <- merge(vars, subst_txt, by = "variable", all.x = TRUE)
vars[, kept_pct := round(100 * kept_in / n_fit)]
vars[is.na(kept_in), `:=`(kept_in = 0L, kept_pct = 0)]
# kept variables first, then dropped; each from most to least often kept
vars[, ord_status := status != "kept"]
setorder(vars, ord_status, -kept_in, full_rank, na.last = TRUE)
vars[, ord_status := NULL]
vars[is.na(instead), instead := ""]
setcolorder(vars, c("variable", "status", "full_rank", "kept_in", "kept_pct", "median_rank", "min_rank",
                    "max_rank", "instead"))
vars[, n_fits := NULL]
setnames(vars, "instead", "when_not_kept_replaced_by")
fwrite(vars, opt$out_variables, sep = "\t", na = "NA", quote = FALSE)

sites <- loo[, {
  s <- variable[kept %in% TRUE]
  .(n_pos_rsq = n_pos_rsq[1], n_kept = length(s),
    same_as_full = setequal(s, full),
    added = paste(sort(setdiff(s, full)), collapse = " "),
    lost = paste(sort(setdiff(full, s)), collapse = " "))
}, by = .(left_out, sitename)]
sites[, ord := suppressWarnings(as.numeric(left_out))]
setorder(sites, -same_as_full, ord, left_out)
sites[, ord := NULL]
sites[, same_as_full := ifelse(same_as_full, "yes", "no")]
fwrite(sites, opt$out_sites, sep = "\t", na = "NA", quote = FALSE)
message(sprintf("%d leave-one-out selections; full-data selection reproduced in %d", n_fit,
                sum(sites$same_as_full == "yes")))
