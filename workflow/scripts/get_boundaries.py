#!/usr/bin/env python
"""Country boundaries clipped to the study area, for map overlays.

--source naturalearth   download Natural Earth 1:10m Admin 0 countries (public domain)
--source <path>         any OGR-readable vector file of boundaries
--source none           write an empty layer (no overlay)
"""
import argparse
import os
import shutil
import sys
import urllib.request

import geopandas as gpd

NE_URL = "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", required=True)
    p.add_argument("--clip", required=True, help="study area GeoJSON")
    p.add_argument("--cache-dir", required=True, help="where the Natural Earth zip is kept")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    area = gpd.read_file(a.clip).to_crs(4326)
    if a.source == "none":
        gpd.GeoDataFrame({"name": []}, geometry=[], crs=4326).to_file(a.out, driver="GeoJSON")
        print("no boundaries requested; wrote empty layer")
        return

    if a.source == "naturalearth":
        os.makedirs(a.cache_dir, exist_ok=True)
        zpath = os.path.join(a.cache_dir, os.path.basename(NE_URL))
        if not os.path.exists(zpath):
            print(f"downloading {NE_URL}")
            req = urllib.request.Request(NE_URL, headers={"User-Agent": "GradientForestMultiSpecies"})
            with urllib.request.urlopen(req, timeout=120) as r, open(zpath + ".part", "wb") as fh:
                shutil.copyfileobj(r, fh)
            os.replace(zpath + ".part", zpath)
        print(f"reading {zpath} ({os.path.getsize(zpath) / 1e6:.1f} MB)")
        gdf = gpd.read_file(f"zip://{zpath}")
        name_col = "NAME" if "NAME" in gdf.columns else gdf.columns[0]
    else:
        gdf = gpd.read_file(a.source)
        name_col = next((c for c in ("NAME", "name", "NAME_0", "ADMIN") if c in gdf.columns), None)

    if gdf.crs is None:
        print("WARNING: boundaries have no CRS; assuming EPSG:4326", file=sys.stderr)
        gdf = gdf.set_crs(4326)
    gdf = gdf.to_crs(4326)
    gdf = gpd.GeoDataFrame({"name": gdf[name_col] if name_col else range(len(gdf))},
                           geometry=gdf.geometry, crs=4326)
    # slightly enlarge the clip so borders run to the map edge
    clip = area.buffer(0.05)
    out = gpd.clip(gdf, clip)
    out = out[~out.geometry.is_empty]
    out.to_file(a.out, driver="GeoJSON")
    print(f"{len(out)} boundary features within the study area: {', '.join(map(str, sorted(out['name'])))}")


if __name__ == "__main__":
    main()
