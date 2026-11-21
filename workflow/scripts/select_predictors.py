#!/usr/bin/env python
"""Choose gradient forest predictors: drop constant variables, then greedily
remove collinear variables (|Pearson r| at the sampled populations above
--threshold). For each offending pair the variable with the higher mean |r|
to the remaining set is dropped, unless it is in --keep-priority."""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--env-sites", required=True)
    p.add_argument("--layers", required=True)
    p.add_argument("--include", default="all", help="'all' or comma-separated variable names")
    p.add_argument("--threshold", default="0.8", help="number, or 'none' to skip pruning")
    p.add_argument("--keep-priority", default="none")
    p.add_argument("--out-list", required=True)
    p.add_argument("--out-report", required=True)
    p.add_argument("--out-heatmap", required=True)
    a = p.parse_args()

    layers = [l.strip() for l in open(a.layers) if l.strip()]
    env = pd.read_csv(a.env_sites, sep="\t")
    vars_ = layers if a.include == "all" else a.include.split(",")
    unknown = set(vars_) - set(layers)
    if unknown:
        raise SystemExit(f"ERROR: predictors not among layers: {sorted(unknown)}")
    keep = set() if a.keep_priority == "none" else set(a.keep_priority.split(","))

    report = []
    X = env[vars_].astype(float)
    for v in vars_:
        if X[v].std() == 0:
            report.append((v, "dropped", "constant across sites", ""))
    current = [v for v in vars_ if X[v].std() > 0]

    if a.threshold != "none":
        thr = float(a.threshold)
        while len(current) > 1:
            r = X[current].corr().abs().to_numpy(copy=True)
            np.fill_diagonal(r, 0)
            i, j = np.unravel_index(np.argmax(r), r.shape)
            if r[i, j] <= thr:
                break
            vi, vj = current[i], current[j]
            mean_r = dict(zip(current, r.sum(axis=0) / (len(current) - 1)))
            cands = [v for v in (vi, vj) if v not in keep] or [vi, vj]
            drop = max(cands, key=lambda v: mean_r[v])
            other = vj if drop == vi else vi
            report.append((drop, "dropped", f"|r|={r[i, j]:.3f} with {other}", f"{mean_r[drop]:.3f}"))
            current.remove(drop)
    for v in current:
        report.append((v, "kept", "", ""))

    rep = pd.DataFrame(report, columns=["variable", "status", "reason", "mean_abs_r"])
    print(f"{len(env)} sites; kept {len(current)} of {len(vars_)} predictors")
    print(rep.to_string(index=False))
    rep.to_csv(a.out_report, sep="\t", index=False)
    with open(a.out_list, "w") as fh:
        fh.write("\n".join(current) + "\n")

    c = X[vars_].corr()
    fig, ax = plt.subplots(figsize=(0.45 * len(vars_) + 2, 0.45 * len(vars_) + 1.5))
    im = ax.imshow(c.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(vars_)), vars_, rotation=90)
    ax.set_yticks(range(len(vars_)), vars_)
    for k, lab in enumerate(ax.get_xticklabels() + ax.get_yticklabels()):
        lab.set_fontweight("bold" if lab.get_text() in current else "normal")
        lab.set_color("black" if lab.get_text() in current else "grey")
    ax.set_title(f"Predictor correlations at {len(env)} sites (bold = kept)")
    fig.colorbar(im, ax=ax, shrink=0.7, label="Pearson r")
    fig.tight_layout()
    fig.savefig(a.out_heatmap, dpi=150)


if __name__ == "__main__":
    main()
