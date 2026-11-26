#!/usr/bin/env python
"""Build a clipped, multiband environmental GeoTIFF from a directory of single-layer
rasters (ESRI ASCII grids, GeoTIFFs, ...).

Layer names are derived from file names with --name-regex / --name-template
(e.g. PIREworldclim-bio_12.asc -> bio12). All layers must share one grid.
Cells that are NA in any layer are set NA in all layers, so every non-NA cell
has a complete predictor vector.
"""
import argparse
import glob
import os
import re
import sys
import warnings

import geopandas as gpd
import numpy as np
import rasterio
import rasterio.mask
from rasterio.crs import CRS


def layer_name(path, regex, template):
    base = os.path.basename(path)
    m = re.search(regex, base, flags=re.IGNORECASE)
    if not m:
        sys.exit(f"ERROR: name regex '{regex}' does not match file {base}")
    groups = m.groups() or (m.group(0),)
    key = groups[0]
    name = template.format(int(key) if key.isdigit() else key)
    return name, (int(key) if key.isdigit() else key)


def block_mean(arr, f):
    """Aggregate a 2-D array by factor f using the mean of non-NaN cells."""
    ny, nx = (arr.shape[0] // f) * f, (arr.shape[1] // f) * f
    a = arr[:ny, :nx].reshape(ny // f, f, nx // f, f)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN blocks -> NaN
        return np.nanmean(a, axis=(1, 3))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--env-dir", required=True)
    p.add_argument("--pattern", default="*.asc")
    p.add_argument("--name-regex", default=r"bio_?(\d+)")
    p.add_argument("--name-template", default="bio{0}")
    p.add_argument("--clip", required=True, help="study area vector (GeoJSON)")
    p.add_argument("--aggregate", type=int, default=1)
    p.add_argument("--out-tif", required=True)
    p.add_argument("--out-layers", required=True)
    a = p.parse_args()

    files = sorted(glob.glob(os.path.join(a.env_dir, a.pattern)))
    files = [f for f in files if not f.endswith((".aux.xml", ".prj"))]
    if not files:
        sys.exit(f"ERROR: no files matching {a.pattern} in {a.env_dir}")
    named = [(layer_name(f, a.name_regex, a.name_template), f) for f in files]
    names = [n for (n, _), _ in named]
    dups = {n for n in names if names.count(n) > 1}
    if dups:
        sys.exit(f"ERROR: duplicate layer names after renaming: {sorted(dups)}")
    named.sort(key=lambda x: (isinstance(x[0][1], str), x[0][1]))

    area = gpd.read_file(a.clip)
    ref = None
    layers, arrays = [], []
    for (name, _), path in named:
        with rasterio.open(path) as src:
            crs = src.crs or CRS.from_epsg(4326)
            if src.crs is None:
                print(f"WARNING: {path} has no CRS; assuming EPSG:4326", file=sys.stderr)
            grid = (src.width, src.height, tuple(round(v, 9) for v in src.transform[:6]))
            if ref is None:
                ref = (grid, crs)
                shapes = list(area.to_crs(crs).geometry)
            elif grid != ref[0] or crs != ref[1]:
                sys.exit(f"ERROR: {path} is on a different grid/CRS from the first layer")
            data, transform = rasterio.mask.mask(src, shapes, crop=True, filled=False)
            arr = data[0].astype("float32").filled(np.nan)
            if src.nodata is not None:
                arr[np.isclose(arr, src.nodata)] = np.nan
            arr[~np.isfinite(arr) | (np.abs(arr) > 1e30)] = np.nan
        layers.append(name)
        arrays.append(arr)
        print(f"{name}\t{os.path.basename(path)}\tmin={np.nanmin(arr):.4g}\tmax={np.nanmax(arr):.4g}")

    stack = np.stack(arrays)
    if a.aggregate > 1:
        stack = np.stack([block_mean(x, a.aggregate) for x in stack])
        transform = transform @ rasterio.Affine.scale(a.aggregate)
    nodata_any = np.isnan(stack).any(axis=0)
    stack[:, nodata_any] = np.nan
    print(f"grid {stack.shape[2]} x {stack.shape[1]}, {int((~nodata_any).sum())} complete cells, "
          f"{len(layers)} layers")

    profile = dict(driver="GTiff", width=stack.shape[2], height=stack.shape[1], count=len(layers),
                   dtype="float32", crs=ref[1], transform=transform, nodata=np.nan,
                   compress="deflate", tiled=True, blockxsize=256, blockysize=256)
    with rasterio.open(a.out_tif, "w", **profile) as dst:
        dst.write(stack)
        for i, n in enumerate(layers, start=1):
            dst.set_band_description(i, n)
    with open(a.out_layers, "w") as fh:
        fh.write("\n".join(layers) + "\n")


if __name__ == "__main__":
    main()
