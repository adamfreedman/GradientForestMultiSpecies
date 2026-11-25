#!/usr/bin/env Rscript
# Per-SNP empirical p-values: how often a randomized fit (environment rows
# permuted) gives a SNP an R^2 at least as large as in the observed model.
#   p = (1 + no. randomized fits with R^2 >= observed) / (n_random + 1)
# gradientForest only reports SNPs with R^2 > 0, so a SNP absent from a fit has
# R^2 = 0 there. Only SNPs in the observed model (R^2 > 0) are tested: the null
# is that their environmental association is no stronger than in randomized
# data. p_bh is the Benjamini-Hochberg adjustment across those SNPs. SNPs with
# observed R^2 <= 0 are kept in the table with NA p-values.
# CHROM, POS and VCF_ID come from genotypes/snp_info.tsv.
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))

opt <- parse_args()
need(opt, c("rsq", "popfreq", "snp_info", "null", "out"))
snps <- setdiff(names(fread(opt$popfreq, nrows = 0)), "popid")
n_rand <- nrow(fread(opt$null, select = "rep"))

res <- data.table(snp = snps, rsq = 0)
obs <- fread(opt$rsq)
if (length(setdiff(obs$snp, snps))) stop("snp_rsq.tsv has SNPs that are not in popfreq.tsv")
res[obs, rsq := i.rsq, on = "snp"]

# decompress with gzip (fread's own .gz support needs R.utils, which is not in the env)
read_tsv <- function(f) if (grepl("\\.gz$", f)) fread(cmd = paste("gzip -dc", shQuote(f))) else fread(f)
nl <- rbindlist(lapply(opt$positional, read_tsv))
if (uniqueN(nl$rep) > n_rand) stop("more replicates in the per-SNP files than in random_null.tsv")
# null R^2 >= observed, for SNPs with observed R^2 > 0 (absent = 0 < observed)
m <- res[rsq > 0][nl, on = "snp", nomatch = NULL]
cnt <- m[i.rsq >= rsq, .N, by = snp]
res[, n_null_ge := 0L]
res[cnt, n_null_ge := i.N, on = "snp"]
res[rsq <= 0, n_null_ge := NA_integer_]
res[, n_random := n_rand]
res[, p_empirical := (1 + n_null_ge) / (1 + n_rand)]
res[, p_bh := p.adjust(p_empirical, method = "BH")]  # NAs are left out of the adjustment
info <- fread(opt$snp_info, colClasses = list(character = c("snp", "CHROM", "VCF_ID")))
if (!all(res$snp %in% info$snp)) stop("snp_info.tsv is missing SNPs that are in popfreq.tsv")
res <- info[res, on = "snp"]
setcolorder(res, c("snp", "CHROM", "POS", "VCF_ID"))
setorderv(res, c("p_empirical", "rsq"), order = c(1L, -1L), na.last = TRUE)
fwrite(res, opt$out, sep = "\t", na = "NA")

message(sprintf("%d SNPs in popfreq, %d tested (observed R^2 > 0), %d randomized fits; smallest possible p = %.3g",
                nrow(res), sum(res$rsq > 0), n_rand, 1 / (n_rand + 1)))
for (a in c(0.01, 0.05, 0.1))
  message(sprintf("alpha %.2f: %d SNPs with p <= alpha, %d with BH-adjusted p <= alpha",
                  a, sum(res$p_empirical <= a, na.rm = TRUE), sum(res$p_bh <= a, na.rm = TRUE)))
