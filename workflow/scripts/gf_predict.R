#!/usr/bin/env Rscript
# Transform the environmental grid with the fitted gradient forest (genomic
# turnover space), then PCA -> RGB colours as in Fitzpatrick & Keller / skinkforest.R.
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(this_dir, "gf_common.R"))
suppressPackageStartupMessages({
  library(gradientForest)
  library(terra)
})

opt <- parse_args()
need(opt, c("model", "env_tif", "env_sites", "pca_vars", "out_turnover", "out_rgb", "out_pca",
            "out_grid", "out_site_pcs", "out_loadings"))

gf <- readRDS(opt$model)
preds <- colnames(gf$X)
env <- rast(opt$env_tif)
missing <- setdiff(preds, names(env))
if (length(missing)) stop("predictors not in env stack: ", paste(missing, collapse = ", "))

grid <- as.data.frame(env[[preds]], xy = TRUE, cells = TRUE, na.rm = TRUE)
message(nrow(grid), " grid cells")
trans <- predict(gf, grid[, preds, drop = FALSE])   # extrapolates linearly beyond sampled range

# turnover rasters: one band per predictor, same grid as the env stack
tmpl <- env[[1]]
to_rast <- function(v, name) {
  r <- rast(tmpl)
  r[grid$cell] <- v
  names(r) <- name
  r
}
turn <- rast(lapply(preds, function(p) to_rast(trans[[p]], p)))
writeRaster(turn, opt$out_turnover, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))

# PCA on the transformed predictors (unscaled: transformed units are R^2-weighted already)
imp <- importance(gf)
pca_vars <- if (identical(opt$pca_vars, "all")) names(imp) else names(imp)[seq_len(min(length(imp), as.integer(opt$pca_vars)))]
if (length(pca_vars) < 2) stop("need at least 2 predictors for the PCA")
pcs <- prcomp(trans[, pca_vars, drop = FALSE])
# sign convention: the largest-magnitude loading of each PC is positive
for (j in seq_len(ncol(pcs$rotation))) {
  if (pcs$rotation[which.max(abs(pcs$rotation[, j])), j] < 0) {
    pcs$rotation[, j] <- -pcs$rotation[, j]
    pcs$x[, j] <- -pcs$x[, j]
  }
}
x <- pcs$x
pc3 <- if (ncol(x) >= 3) x[, 3] else 0
rgb_raw <- cbind(r = x[, 1] + x[, 2], g = -x[, 2], b = pc3 + x[, 2] - x[, 1])
rng <- apply(rgb_raw, 2, range)
scale255 <- function(v, r) if (diff(r) > 0) (v - r[1]) / diff(r) * 255 else rep(127, length(v))
rgb_val <- sapply(1:3, function(i) scale255(rgb_raw[, i], rng[, i]))
colnames(rgb_val) <- c("r", "g", "b")

# 0 is reserved for NA in the 8-bit RGB GeoTIFF
rgbr <- rast(lapply(1:3, function(i) to_rast(pmax(1, round(rgb_val[, i])), c("red", "green", "blue")[i])))
writeRaster(rgbr, opt$out_rgb, overwrite = TRUE, datatype = "INT1U", NAflag = 0,
            gdal = c("COMPRESS=DEFLATE", "PHOTOMETRIC=RGB"))

npc <- min(3, ncol(x))
grid_out <- data.frame(x = grid$x, y = grid$y, x[, seq_len(npc), drop = FALSE], rgb_val)
saveRDS(grid_out, opt$out_grid)

# populations in the same PC space (transformed site values from the training data)
site_trans <- predict(gf)
site_pc <- predict(pcs, site_trans[, pca_vars, drop = FALSE])
sites <- fread(opt$env_sites, colClasses = list(character = "popid"))
sites <- sites[match(attr(gf, "popid"), sites$popid)]  # gf$X row order
site_rgb <- sapply(1:3, function(i) {
  v <- cbind(site_pc[, 1] + site_pc[, 2], -site_pc[, 2],
             (if (ncol(site_pc) >= 3) site_pc[, 3] else 0) + site_pc[, 2] - site_pc[, 1])[, i]
  pmin(255, pmax(0, scale255(v, rng[, i])))
})
site_out <- data.table(popid = sites$popid, sitename = sites$sitename, lat = sites$lat, long = sites$long,
                       site_pc[, seq_len(npc), drop = FALSE], r = site_rgb[, 1], g = site_rgb[, 2], b = site_rgb[, 3])
fwrite(site_out, opt$out_site_pcs, sep = "\t")

load <- data.table(variable = rownames(pcs$rotation), pcs$rotation[, seq_len(npc), drop = FALSE])
fwrite(load, opt$out_loadings, sep = "\t")
ve <- pcs$sdev^2 / sum(pcs$sdev^2)
saveRDS(list(pca = pcs[c("sdev", "rotation", "center", "scale")], pca_vars = pca_vars,
             var_explained = ve, rgb_range = rng, importance = imp), opt$out_pca)
message("PC variance explained: ", paste(sprintf("%.1f%%", 100 * ve[seq_len(npc)]), collapse = ", "))
