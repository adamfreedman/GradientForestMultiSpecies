"""Shared helpers for the Python workflow scripts."""
import re
import sys

import numpy as np
import pandas as pd
import rasterio

KM_PER_DEG = 111.32


def read_popmap(path):
    """Two-column table (sampleid popid) with a header. Comma delimited if the header
    contains a comma, otherwise any run of spaces and/or tabs separates columns."""
    with open(path) as fh:
        header = fh.readline()
    sep = "," if "," in header else r"\s+"
    # index_col=False: a trailing delimiter must not turn the first column into the index
    df = pd.read_csv(path, sep=sep, engine="python", dtype=str, index_col=False)
    df.columns = [c.strip().lower() for c in df.columns]
    if not {"sampleid", "popid"} <= set(df.columns):
        sys.exit(f"ERROR: {path} must have columns 'sampleid' and 'popid' (found {list(df.columns)})")
    df = df[["sampleid", "popid"]].apply(lambda s: s.str.strip())
    empty = df.isna().any(axis=1) | (df == "").any(axis=1)
    if empty.any():
        rows = (df.index[empty] + 2).tolist()  # file line numbers (header is line 1)
        sys.exit(f"ERROR: {path}: empty sampleid or popid on line(s) {rows[:10]}"
                 f"{' ...' if len(rows) > 10 else ''}")
    return df


def read_latlongs(path):
    """CSV with popid, sitename, lat, long (long/lon/longitude accepted)."""
    df = pd.read_csv(path, sep=None, engine="python", dtype={"popid": str})
    df.columns = [c.strip().lower() for c in df.columns]
    rename = {"lon": "long", "longitude": "long", "latitude": "lat", "site": "sitename"}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    missing = {"popid", "lat", "long"} - set(df.columns)
    if missing:
        sys.exit(f"ERROR: {path} is missing columns {sorted(missing)}")
    if "sitename" not in df.columns:
        df["sitename"] = df["popid"]
    df["popid"] = df["popid"].str.strip()
    return df[["popid", "sitename", "lat", "long"]]


def normalize_ids(ids, regex):
    return [re.sub(regex, "", i) for i in ids] if regex else list(ids)


class EnvStack:
    """Multiband environmental raster held in memory (small enough at regional scale)."""

    def __init__(self, path):
        with rasterio.open(path) as src:
            self.data = src.read()  # (bands, rows, cols)
            self.transform = src.transform
            self.bounds = src.bounds
            self.names = [d or f"band{i + 1}" for i, d in enumerate(src.descriptions)]
        self.valid = ~np.isnan(self.data).any(axis=0)
        self._vxy = None

    def rowcol(self, lon, lat):
        col, row = ~self.transform * (lon, lat)
        return int(np.floor(row)), int(np.floor(col))

    def inside(self, lon, lat):
        b = self.bounds
        return b.left <= lon < b.right and b.bottom < lat <= b.top

    def valid_xy(self):
        if self._vxy is None:
            rows, cols = np.nonzero(self.valid)
            xs, ys = self.transform * (cols + 0.5, rows + 0.5)
            self._vxy = (rows, cols, np.asarray(xs), np.asarray(ys))
        return self._vxy

    def nearest_valid(self, lon, lat):
        """Nearest complete cell; returns (row, col, distance_km)."""
        rows, cols, xs, ys = self.valid_xy()
        d = haversine_km(lon, lat, xs, ys)
        i = int(np.argmin(d))
        return rows[i], cols[i], float(d[i])

    def sample(self, lon, lat, buffer_km=0.0, max_snap_km=0.0):
        """Return (values, snapped_km, n_cells). values is NaN if nothing usable."""
        nan = np.full(len(self.names), np.nan)
        if not self.inside(lon, lat):
            return nan, np.nan, 0
        if buffer_km and buffer_km > 0:
            rows, cols, xs, ys = self.valid_xy()
            sel = haversine_km(lon, lat, xs, ys) <= buffer_km
            if sel.any():
                return self.data[:, rows[sel], cols[sel]].mean(axis=1), 0.0, int(sel.sum())
        r, c = self.rowcol(lon, lat)
        if self.valid[r, c]:
            return self.data[:, r, c], 0.0, 1
        if max_snap_km and max_snap_km > 0:
            r, c, d = self.nearest_valid(lon, lat)
            if d <= max_snap_km:
                return self.data[:, r, c], d, 1
            return nan, d, 0
        return nan, np.nan, 0


def haversine_km(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, (lon1, lat1, lon2, lat2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 6371.0088 * 2 * np.arcsin(np.sqrt(a))
