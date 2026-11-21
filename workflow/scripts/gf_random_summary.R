#!/usr/bin/env Rscript
# Gather randomization replicates; empirical p-values for the observed statistics.
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))

opt <- parse_args()
need(opt, c("summary", "importance", "out_null", "out_pvals"))
null <- rbindlist(lapply(opt$positional, fread), fill = TRUE)
setorder(null, rep)
fwrite(null, opt$out_null, sep = "\t")

obs <- fread(opt$summary)
imp <- fread(opt$importance)
emp_p <- function(o, n) (1 + sum(n >= o, na.rm = TRUE)) / (1 + sum(!is.na(n)))

res <- rbind(
  data.table(statistic = "n_pos_rsq", observed = obs$n_pos_rsq),
  data.table(statistic = "mean_rsq", observed = obs$mean_rsq),
  data.table(statistic = "sum_rsq", observed = obs$sum_rsq),
  data.table(statistic = paste0("imp_", imp$variable), observed = imp$importance),
  data.table(statistic = paste0("acc_", imp$variable), observed = imp$accuracy)
)
res[, null_mean := sapply(statistic, function(s) mean(null[[s]], na.rm = TRUE))]
res[, null_q95 := sapply(statistic, function(s) quantile(null[[s]], 0.95, na.rm = TRUE))]
res[, p_empirical := mapply(function(s, o) emp_p(o, null[[s]]), statistic, observed)]
res[, n_random := nrow(null)]
fwrite(res, opt$out_pvals, sep = "\t")
print(res)
