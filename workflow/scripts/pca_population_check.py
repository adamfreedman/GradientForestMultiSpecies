#!/usr/bin/env python
"""Flag individuals whose genotypes place them with another population (genotype PCA QC).

Uses each individual's coordinates on the first --n-pcs PCs (sample_scores.tsv).
For each individual, every population's median is computed without that
individual (so it cannot pull its own population's median towards itself), and:
  dist_own            distance to its own population's median
  nearest_other_pop   the other population with the closest median
  dist_nearest_other  distance to that median
  ratio               dist_own / dist_nearest_other (> 1: closer to another population)
  spread_ratio        dist_own / the median distance of the other members of its
                      population to their median (how far outside its own
                      population's usual spread it lies)
flagged = yes when ratio >= --min-ratio AND spread_ratio >= --min-spread-ratio.
Both are needed: populations that are genetically similar overlap in the PCA, so
many individuals sit slightly closer to a neighbouring population's median without
being out of place. When an individual's population has fewer than 2 other
members, its spread cannot be estimated and the median spread over all populations
with at least 3 individuals is used instead (note column); an individual alone in
its population is not checked.
A flagged individual may be mislabelled in the popmap, a sample mix-up or
contamination, a migrant or admixed individual, or a population that is
genetically heterogeneous.
"""
import argparse

import numpy as np
import pandas as pd


def check(scores, n_pcs=4, min_ratio=2.0, min_spread_ratio=3.0):
    pcs = [c for c in scores.columns if c.startswith("PC")][:n_pcs]
    X = scores[pcs].to_numpy(float)
    pops = scores.popid.to_numpy()
    # typical spread of a population, for populations too small to estimate their own:
    # median, over populations with >= 3 individuals, of the members' median distance
    # to their population median
    own_spread = []
    for p in np.unique(pops):
        m = X[pops == p]
        if len(m) >= 3:
            own_spread.append(np.median(np.linalg.norm(m - np.median(m, axis=0), axis=1)))
    pooled_spread = float(np.median(own_spread)) if own_spread else np.nan
    rows = []
    for i in range(len(scores)):
        others = np.arange(len(scores)) != i
        med = {p: np.median(X[(pops == p) & others], axis=0) for p in np.unique(pops[others])}
        own = pops[i]
        row = dict(sample_id=scores.sample_id.iloc[i], popid=own, nearest_other_pop="", dist_own=np.nan,
                   dist_nearest_other=np.nan, ratio=np.nan, spread_ratio=np.nan, flagged="no", note="")
        if own not in med:
            row["note"] = "only individual in its population: not checked"
            rows.append(row)
            continue
        dist = {p: float(np.linalg.norm(X[i] - m)) for p, m in med.items()}
        other = {p: d for p, d in dist.items() if p != own}
        if not other:
            row["note"] = "only one population: not checked"
            rows.append(row)
            continue
        near = min(other, key=other.get)
        mates = X[(pops == own) & others]
        if len(mates) >= 2:
            spread = float(np.median(np.linalg.norm(mates - med[own], axis=1)))
        else:
            spread = pooled_spread
            row["note"] = "fewer than 3 individuals in its population: spread of all populations used"
        row.update(nearest_other_pop=near, dist_own=dist[own], dist_nearest_other=other[near],
                   ratio=dist[own] / other[near] if other[near] > 0 else np.inf,
                   spread_ratio=dist[own] / spread if spread and spread > 0 else np.nan)
        spread_ok = np.isnan(row["spread_ratio"]) or row["spread_ratio"] >= min_spread_ratio
        if row["ratio"] >= min_ratio and spread_ok:
            row["flagged"] = "yes"
        rows.append(row)
    return pd.DataFrame(rows), len(pcs)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scores", required=True, help="pca/sample_scores.tsv")
    p.add_argument("--n-pcs", type=int, default=4)
    p.add_argument("--min-ratio", type=float, default=2.0)
    p.add_argument("--min-spread-ratio", type=float, default=3.0)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    scores = pd.read_csv(a.scores, sep="\t", dtype={"sample_id": str, "popid": str})
    res, k = check(scores, a.n_pcs, a.min_ratio, a.min_spread_ratio)
    res.to_csv(a.out, sep="\t", index=False, float_format="%.4g")
    flagged = res[res.flagged == "yes"]
    print(f"PCs used: {k}; flagged {len(flagged)} of {len(res)} individuals "
          f"(ratio >= {a.min_ratio:g} and spread_ratio >= {a.min_spread_ratio:g})")
    if len(flagged):
        print(flagged.to_string(index=False))


if __name__ == "__main__":
    main()
