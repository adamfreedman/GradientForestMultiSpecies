#!/usr/bin/env python
"""Environmental values at each population location (cell value, or the mean of
cells within --buffer-km; NA cells snap to the nearest valid cell within
--max-snap-km)."""
import argparse
import sys

import numpy as np
import pandas as pd

from gfutils import EnvStack, read_latlongs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--latlongs", required=True)
    p.add_argument("--sample-pop", required=True)
    p.add_argument("--env-tif", required=True)
    p.add_argument("--buffer-km", type=float, default=0)
    p.add_argument("--max-snap-km", type=float, default=0)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    sp = pd.read_csv(a.sample_pop, sep="\t", dtype=str)
    counts = sp.popid.value_counts()
    ll = read_latlongs(a.latlongs)
    ll = ll[ll.popid.isin(counts.index)].copy()
    env = EnvStack(a.env_tif)

    rows = []
    for _, r in ll.iterrows():
        vals, snap, ncell = env.sample(float(r["long"]), float(r["lat"]), a.buffer_km, a.max_snap_km)
        rows.append([r.popid, r.sitename, r.lat, r["long"], int(counts[r.popid]), snap, ncell, *vals])
    out = pd.DataFrame(rows, columns=["popid", "sitename", "lat", "long", "n_ind", "snapped_km", "n_cells",
                                      *env.names])
    bad = out[out[env.names].isna().any(axis=1)]
    print(out.to_string(index=False))
    if len(bad):
        sys.exit(f"ERROR: no environmental data for populations {list(bad.popid)}")
    out.to_csv(a.out, sep="\t", index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
