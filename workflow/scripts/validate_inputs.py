#!/usr/bin/env python
"""Cross-check a species' VCF, population map and population coordinates.

Individuals that do not map to a population with coordinates are excluded
with a warning and listed in --out-unmapped:
  * in the VCF but not in the popmap
  * in the popmap but not in the VCF
  * in both, but their popid has no row in latlongs (e.g. a population left
    out because it has too few individuals)

Fails (exit 1, errors on stderr, no outputs) on: a latlongs population with no
individuals in the VCF, duplicated IDs, populations outside the environmental
grid or on NA cells farther than --max-snap-km from valid data.
Writes a report, the vcf_id -> popid table and the unmapped-individuals table.
"""
import argparse
import os
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
    p.add_argument("--out-report", required=True)
    p.add_argument("--out-sample-pop", required=True)
    p.add_argument("--out-unmapped", required=True)
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

    ll = read_latlongs(a.latlongs)
    info.append(f"latlongs populations: {len(ll)}")
    if ll.popid.duplicated().any():
        errors.append(f"duplicate popids in latlongs: {sorted(ll.popid[ll.popid.duplicated()])}")
    ll_pops = set(ll.popid)
    vcf_set = set(norm)

    # individuals that do not map to a population with coordinates
    unmapped = []  # sample_id, popid, missing_in_vcf, missing_in_popmap, popid_in_latlongs, reason
    for v, n in zip(vcf_ids, norm):
        if n not in pm_d:
            unmapped.append((v, "", "no", "yes", "", "in VCF, not in popmap"))
        elif pm_d[n] not in ll_pops:
            unmapped.append((v, pm_d[n], "no", "no", "no", "population has no coordinates in latlongs"))
    for sid in pm.sampleid.drop_duplicates():
        if sid not in vcf_set:
            pop = pm_d[sid]
            unmapped.append((sid, pop, "yes", "no", "yes" if pop in ll_pops else "no", "in popmap, not in VCF"))

    def short(ids, k=10):
        return ", ".join(ids[:k]) + (f", ... ({len(ids) - k} more)" if len(ids) > k else "")

    by_reason = {}
    for r in unmapped:
        by_reason.setdefault(r[5], []).append(r)
    if "in VCF, not in popmap" in by_reason:
        ids = [r[0] for r in by_reason["in VCF, not in popmap"]]
        warns.append(f"{len(ids)} VCF samples not in popmap (excluded): {short(ids)}")
    if "population has no coordinates in latlongs" in by_reason:
        rs = by_reason["population has no coordinates in latlongs"]
        pops = sorted({r[1] for r in rs}, key=lambda x: (len(x), x))
        counts = ", ".join(f"{pop} (n={sum(r[1] == pop for r in rs)})" for pop in pops)
        warns.append(f"{len(rs)} VCF samples in popids with no coordinates in latlongs (excluded): {counts}")
    if "in popmap, not in VCF" in by_reason:
        ids = [r[0] for r in by_reason["in popmap, not in VCF"]]
        warns.append(f"{len(ids)} popmap samples not in VCF (ignored): {short(ids)}")
    if unmapped:
        warns.append(f"{len(unmapped)} individuals in total do not map to a population; listed in "
                     f"{os.path.basename(a.out_unmapped)}")

    mapped = [(v, pm_d[n]) for v, n in zip(vcf_ids, norm) if n in pm_d and pm_d[n] in ll_pops]
    info.append(f"VCF samples mapped to a population with coordinates: {len(mapped)}")
    pops_used = sorted({pop for _, pop in mapped}, key=lambda x: (len(x), x))
    empty = ll[~ll.popid.isin(pops_used)]
    for _, row in empty.iterrows():
        n_pm = int((pm.popid == row.popid).sum())
        errors.append(f"latlongs population {row.popid} ({row.sitename}) has no individuals in the VCF "
                      f"({n_pm} in popmap); remove it from latlongs or fix the popmap/VCF")

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
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(f"validation FAILED with {len(errors)} error(s)")
    with open(a.out_report, "w") as fh:
        fh.write(report)
    with open(a.out_sample_pop, "w") as fh:
        fh.write("vcf_id\tpopid\n")
        for v, pop in mapped:
            fh.write(f"{v}\t{pop}\n")
    with open(a.out_unmapped, "w") as fh:
        fh.write("sample_id\tpopid\tmissing_in_vcf\tmissing_in_popmap\tpopid_in_latlongs\treason\n")
        for r in unmapped:
            fh.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main()
