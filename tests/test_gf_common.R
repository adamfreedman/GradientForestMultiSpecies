#!/usr/bin/env Rscript
# Tests for workflow/scripts/gf_common.R. Needs data.table and gradientForest,
# i.e. the workflow's R conda environment. Run from the repository root:
#   Rscript tests/test_gf_common.R
source("workflow/scripts/gf_common.R")

n_fail <- 0L
check <- function(desc, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { message("  error: ", conditionMessage(e)); FALSE })
  cat(if (ok) "ok  " else "FAIL", desc, "\n")
  if (!ok) n_fail <<- n_fail + 1L
}
errors_with <- function(expr, pattern) {
  msg <- tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))
  grepl(pattern, msg)
}

# argument parsing --------------------------------------------------------------
o <- parse_args(c("--a=1", "--b=x=y", "pos1", "--empty=", "pos2"))
check("parse_args key=value", identical(o$a, "1") && identical(o$b, "x=y") && identical(o$empty, ""))
check("parse_args positional", identical(o$positional, c("pos1", "pos2")))
check("need() reports missing arguments", errors_with(need(o, c("a", "zz")), "--zz"))

# maxLevel ------------------------------------------------------------------------
check("maxLevel auto, 14 sites", resolve_maxlevel("auto", 14) == 1)
check("maxLevel auto, 100 sites", resolve_maxlevel("auto", 100) == 4)
check("maxLevel auto never below 1", resolve_maxlevel("auto", 3) == 1)
check("maxLevel fixed", resolve_maxlevel("3", 100) == 3)

# importance floor ----------------------------------------------------------------
check("floor absent", is.na(parse_floor(list(), "importance")))
check("floor none", is.na(parse_floor(list(importance_floor = "none"), "importance")))
check("floor value", parse_floor(list(importance_floor = "0.01"), "importance") == 0.01)
check("floor >= 1 rejected", errors_with(parse_floor(list(importance_floor = "1"), "importance"), "fraction"))
check("floor not a number", errors_with(parse_floor(list(importance_floor = "abc"), "importance"), "fraction"))
check("floor needs a fit-based ranking",
      errors_with(parse_floor(list(importance_floor = "0.01"), "correlation"), "importance or accuracy"))

# selection by correlation (no model fit) ----------------------------------------
set.seed(1)
n <- 12
X <- data.frame(a = 1:n + rnorm(n, 0, 0.5))
X$b <- X$a + rnorm(n, 0, 0.5)        # |r| with a > 0.8
X$c <- rnorm(n)                      # independent
X$d <- X$c * 0.2 + rnorm(n)          # weakly related to c
r <- abs(cor(X))
check("test data: a and b correlated, others not",
      r["a", "b"] > 0.8 && max(r[c("a", "b"), c("c", "d")]) < 0.8 && r["c", "d"] < 0.8)
s <- select_predictors(X, NULL, list(), 0.8, "correlation")
check("correlation: one of a/b kept, c and d kept",
      sum(c("a", "b") %in% s$selected) == 1 && all(c("c", "d") %in% s$selected))
dropped <- setdiff(c("a", "b"), s$selected)
check("correlation: reason names the kept partner", grepl("^\\|r\\|=", s$reason[dropped]) &&
        grepl(setdiff(c("a", "b"), dropped), s$reason[dropped]))
check("no two kept variables above the threshold",
      { k <- s$selected; m <- r[k, k]; diag(m) <- 0; max(m) <= 0.8 })
check("threshold NA keeps everything", setequal(select_predictors(X, NULL, list(), NA, "correlation")$selected,
                                                names(X)))
s2 <- select_predictors(X, NULL, list(), 0.8, "correlation", keep_priority = dropped)
check("keep_priority wins within a correlated pair", dropped %in% s2$selected && s2$rank[[dropped]] == 1)

# selection by a screening gradient forest ----------------------------------------
# SNP frequencies follow b and, more weakly, c; d is noise. (a and b are too
# correlated for the fit to tell which one drives the SNPs.)
set.seed(2)
Y <- as.data.frame(sapply(1:60, function(i) {
  z <- if (i <= 40) X$b else X$c
  plogis(as.numeric(scale(z)) * runif(1, 2, 3) + rnorm(n, 0, 0.2))
}))
opt <- list(screen_ntree = "200", maxLevel = "auto", corr_threshold = "0.5", nbin = "101")
s3 <- select_predictors(X, Y, opt, 0.8, "importance")
kept_ab <- intersect(c("a", "b"), s3$selected)
check("importance: of the correlated pair, only the higher-scoring one is kept",
      length(kept_ab) == 1 && s3$score[[kept_ab]] == max(s3$score[c("a", "b")]))
check("importance: the noise variable d ranks below the signal variables",
      s3$rank[["d"]] > s3$rank[[kept_ab]] && s3$rank[["d"]] > s3$rank[["c"]])
check("importance: shares of the positive scores sum to 1",
      all(is.finite(s3$score)) && abs(sum(s3$share[s3$score > 0]) - 1) < 1e-9)
check("importance: some SNPs enter the screening model", s3$n_pos_rsq > 0)
fl <- sort(s3$share[s3$selected])[2] + 1e-6   # floor just above the second-smallest kept share
s4 <- select_predictors(X, Y, opt, 0.8, "importance", importance_floor = fl)
check("floor drops low-share variables only", all(s4$share[s4$selected] >= fl) &&
        all(grepl("below importance floor", s4$reason[setdiff(s3$selected, s4$selected)])))
check("floor never removes the top-ranked variable", names(which.max(s3$score)) %in% s4$selected)
low <- names(which.min(s3$share[s3$selected]))
s5 <- select_predictors(X, Y, opt, 0.8, "importance", keep_priority = low, importance_floor = fl)
check("keep_priority variables are exempt from the floor", low %in% s5$selected)

cat(if (n_fail) sprintf("\n%d test(s) FAILED\n", n_fail) else "\nall tests passed\n")
quit(status = if (n_fail) 1L else 0L)
