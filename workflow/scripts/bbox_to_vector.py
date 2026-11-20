#!/usr/bin/env python
"""Convert a study-area polygon (shapefile, GeoPackage, GeoJSON, ...) or a plain
bounding box into GeoJSON (WGS84) plus a one-line xmin/ymin/xmax/ymax table.

Usable on its own, without GIS software:
    python bbox_to_vector.py --vector PIREbb.shp --geojson study_area.geojson --bbox bbox.tsv
    python bbox_to_vector.py --extent 7.7 -1.6 19.1 13.4 --geojson study_area.geojson
"""
import argparse
import sys

import geopandas as gpd
from shapely.geometry import box


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--vector", help="any OGR-readable vector file")
    src.add_argument("--extent", nargs=4, type=float, metavar=("XMIN", "YMIN", "XMAX", "YMAX"),
                     help="bounding box in decimal degrees (WGS84)")
    p.add_argument("--geojson", required=True, help="output GeoJSON")
    p.add_argument("--bbox", help="output TSV with xmin ymin xmax ymax")
    a = p.parse_args()

    if a.vector:
        gdf = gpd.read_file(a.vector)
        if gdf.crs is None:
            print("WARNING: input has no CRS; assuming EPSG:4326", file=sys.stderr)
            gdf = gdf.set_crs(4326)
        print(f"read {len(gdf)} feature(s) from {a.vector}, CRS {gdf.crs.to_string()}")
        gdf = gdf.to_crs(4326)
        # dissolve to a single (multi)polygon so downstream masking is simple
        gdf = gpd.GeoDataFrame({"name": ["study_area"]}, geometry=[gdf.union_all()], crs=4326)
    else:
        xmin, ymin, xmax, ymax = a.extent
        if xmin >= xmax or ymin >= ymax:
            sys.exit("ERROR: extent must be xmin ymin xmax ymax with min < max")
        gdf = gpd.GeoDataFrame({"name": ["study_area"]}, geometry=[box(xmin, ymin, xmax, ymax)], crs=4326)

    if not gdf.geometry.iloc[0].is_valid:
        print("WARNING: geometry invalid; repairing with buffer(0)", file=sys.stderr)
        gdf.geometry = gdf.geometry.buffer(0)

    gdf.to_file(a.geojson, driver="GeoJSON")
    xmin, ymin, xmax, ymax = gdf.total_bounds
    print(f"bounds: xmin={xmin} ymin={ymin} xmax={xmax} ymax={ymax}")
    print(f"area is rectangle: {gdf.geometry.iloc[0].equals(box(xmin, ymin, xmax, ymax))}")
    if a.bbox:
        with open(a.bbox, "w") as fh:
            fh.write("xmin\tymin\txmax\tymax\n")
            fh.write(f"{xmin}\t{ymin}\t{xmax}\t{ymax}\n")


if __name__ == "__main__":
    main()
