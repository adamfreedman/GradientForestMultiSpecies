#!/usr/bin/env python
"""Pick the predictors whose observed importance beats the randomized fits.

A predictor is kept if its empirical p (from gf/random_pvalues.tsv) is <= alpha
for any of the requested metrics (importance = R^2-weighted, accuracy).
"""
import argparse
import sys

import pandas as pd

PREFIX = {"importance": "imp_", "accuracy": "acc_"}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pvals", required=True)
    p.add_argument("--alpha", type=float, required=True)
    p.add_argument("--metrics", required=True, help="comma-separated: importance,accuracy")
    p.add_argument("--min-predictors", type=int, default=2)
    p.add_argument("--out-list", required=True)
    p.add_argument("--out-report", required=True)
    a = p.parse_args()

    metrics = [m.strip() for m in a.metrics.split(",") if m.strip()]
    bad = [m for m in metrics if m not in PREFIX]
    if bad:
        sys.exit(f"unknown metric(s) {bad}; use {list(PREFIX)}")

    pv = pd.read_csv(a.pvals, sep="\t")
    rows = {}
    for m in metrics:
        sub = pv[pv["statistic"].str.startswith(PREFIX[m])]
        for _, r in sub.iterrows():
            v = r["statistic"][len(PREFIX[m]):]
            rows.setdefault(v, {"variable": v})
            rows[v][f"{m}_observed"] = r["observed"]
            rows[v][f"{m}_p"] = r["p_empirical"]
    rep = pd.DataFrame(list(rows.values()))
    pcols = [f"{m}_p" for m in metrics]
    rep["kept"] = (rep[pcols] <= a.alpha).any(axis=1)
    rep["reason"] = rep.apply(
        lambda r: ", ".join(f"{m} p={r[f'{m}_p']:.3g}" for m in metrics if r[f"{m}_p"] <= a.alpha)
        or f"no metric with p <= {a.alpha}",
        axis=1,
    )
    rep = rep.sort_values(["kept", pcols[0]], ascending=[False, True])
    rep.to_csv(a.out_report, sep="\t", index=False)
    kept = rep.loc[rep["kept"], "variable"].tolist()
    print(rep.to_string(index=False))
    if len(kept) < a.min_predictors:
        sys.exit(
            f"only {len(kept)} predictor(s) pass p <= {a.alpha} ({', '.join(kept) or 'none'}); "
            f"need {a.min_predictors} for the refit. Set gf.refit_significant.enabled: false "
            "or relax alpha / min_predictors."
        )
    with open(a.out_list, "w") as fh:
        fh.write("\n".join(kept) + "\n")
    print(f"kept {len(kept)}: {', '.join(kept)}")


if __name__ == "__main__":
    main()
