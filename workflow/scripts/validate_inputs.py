#!/usr/bin/env python
"""Cross-check a species' VCF, population map and population coordinates.

Fails (exit 1, no outputs) on: VCF samples missing from the popmap (unless
--allow-unmapped), duplicated IDs, popids without coordinates, populations
outside the environmental grid or on NA cells farther than --max-snap-km from
valid data. Writes a report and a clean vcf_id -> popid table.
"""
import argparse
import sys

import numpy as np
from cyvcf2 import VCF

from gfutils import EnvStack, normalize_ids, read_latlongs, read_popmap


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True)
    p.add_argument("--popmap", required=True)
    p.add_argument("--latlongs", required=True)
    p.add_argument("--env-tif", required=True)
    p.add_argument("--id-regex", default="")
    p.add_argument("--max-snap-km", type=float, default=0)
    p.add_argument("--allow-unmapped", action="store_true")
    p.add_argument("--out-report", required=True)
    p.add_argument("--out-sample-pop", required=True)
    a = p.parse_args()

    errors, warns, info = [], [], []

    vcf_ids = list(VCF(a.vcf).samples)
    norm = normalize_ids(vcf_ids, a.id_regex)
    info.append(f"VCF samples: {len(vcf_ids)}")
    if len(set(norm)) != len(norm):
        d = sorted({x for x in norm if norm.count(x) > 1})
        errors.append(f"duplicate VCF sample IDs (after id_regex): {d}")

    pm = read_popmap(a.popmap)
    info.append(f"popmap rows: {len(pm)}")
    if pm.sampleid.duplicated().any():
        errors.append(f"duplicate sampleids in popmap: {sorted(pm.sampleid[pm.sampleid.duplicated()].unique())}")
    pm_d = dict(zip(pm.sampleid, pm.popid))

    unmapped = [v for v, n in zip(vcf_ids, norm) if n not in pm_d]
    if unmapped:
        msg = f"{len(unmapped)} VCF samples not in popmap: {unmapped}"
        (warns if a.allow_unmapped else errors).append(msg)
    absent = sorted(set(pm.sampleid) - set(norm))
    if absent:
        warns.append(f"{len(absent)} popmap samples not in VCF (ignored): {absent}")

    ll = read_latlongs(a.latlongs)
    if ll.popid.duplicated().any():
        errors.append(f"duplicate popids in latlongs: {sorted(ll.popid[ll.popid.duplicated()])}")
    mapped = [(v, pm_d[n]) for v, n in zip(vcf_ids, norm) if n in pm_d]
    pops_used = sorted({pop for _, pop in mapped}, key=lambda x: (len(x), x))
    no_coord = [x for x in pops_used if x not in set(ll.popid)]
    if no_coord:
        errors.append(f"popids in popmap without coordinates in latlongs: {no_coord}")
    unused = sorted(set(ll.popid) - set(pops_used))
    if unused:
        warns.append(f"latlongs populations with no VCF samples (ignored): {unused}")

    env = EnvStack(a.env_tif)
    info.append(f"\npopulation  site  lat  long  n_ind  env_cell")
    for _, row in ll[ll.popid.isin(pops_used)].iterrows():
        n = sum(1 for _, pop in mapped if pop == row.popid)
        lon, lat = float(row["long"]), float(row["lat"])
        if not np.isfinite([lon, lat]).all():
            errors.append(f"pop {row.popid}: non-numeric coordinates")
            continue
        if abs(lat) > 90 or abs(lon) > 180:
            errors.append(f"pop {row.popid}: coordinates out of range (lat={lat}, long={lon}); lat/long swapped?")
            continue
        if not env.inside(lon, lat):
            errors.append(f"pop {row.popid} ({row.sitename}) lies outside the environmental grid / study area")
            status = "OUTSIDE"
        else:
            r, c = env.rowcol(lon, lat)
            if env.valid[r, c]:
                status = "ok"
            else:
                _, _, d = env.nearest_valid(lon, lat)
                status = f"NA cell; nearest valid {d:.2f} km"
                if d <= a.max_snap_km:
                    warns.append(f"pop {row.popid} ({row.sitename}) on NA cell; will snap {d:.2f} km")
                else:
                    errors.append(f"pop {row.popid} ({row.sitename}) on NA cell, nearest valid cell {d:.2f} km "
                                  f"> max_snap_km {a.max_snap_km}")
        info.append(f"{row.popid}  {row.sitename}  {lat}  {lon}  {n}  {status}")

    lines = ["# input validation", *info, "", f"# warnings ({len(warns)})", *warns, "",
             f"# errors ({len(errors)})", *errors]
    report = "\n".join(lines) + "\n"
    print(report)
    if errors:
        sys.exit("validation FAILED; see errors above")
    with open(a.out_report, "w") as fh:
        fh.write(report)
    with open(a.out_sample_pop, "w") as fh:
        fh.write("vcf_id\tpopid\n")
        for v, pop in mapped:
            fh.write(f"{v}\t{pop}\n")


if __name__ == "__main__":
    main()
