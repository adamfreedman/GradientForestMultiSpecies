#!/usr/bin/env python
"""Summarize the ADMIXTURE runs: CV error per K, pong's modes, the choice of K, and
ancestry tables whose cluster columns are aligned across K.

Choice of K: the K with the lowest mean CV error over the replicate runs (K_cv).
Stability: pong groups the runs at each K into modes (runs whose clusters match with
similarity above pong's threshold). K is stable when pong's largest (major) mode holds
at least --stable-min-frac of the runs.

Outputs (in --outdir)
  cv_error.tsv           one row per run: K, rep, seed, cv_error, loglik
  k_summary.tsv          one row per K: CV error (mean, sd, min), pong modes, stable
  selected_k.tsv         one row: K_cv and its stability, and the same for K_cv - 1
  ancestry_K<k>.tsv      per K: one row per individual (sample_id, popid, then
                         cluster1..clusterK), from the representative run of pong's major
                         mode. Clusters are aligned across K: cluster j at K is the
                         cluster that best matches cluster j at K - 1, and clusterK is
                         the one that is new at this K.
  selected_k_modes.tsv   at K_cv: the representative run of every pong mode (largest
                         first), clusters aligned with ancestry_K<K_cv>.tsv
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

CV_RE = re.compile(r"^CV error \(K=(\d+)\):\s*([0-9.eE+-]+)")
LL_RE = re.compile(r"^Loglikelihood:\s*([0-9.eE+-]+)")


def parse_log(path):
    cv = ll = None
    with open(path) as fh:
        for line in fh:
            m = CV_RE.match(line)
            if m:
                cv = float(m.group(2))
            m = LL_RE.match(line)
            if m:
                ll = float(m.group(1))
    if cv is None:
        sys.exit(f"ERROR: no CV error in {path} (did ADMIXTURE finish?)")
    return cv, ll


def parse_pong(path):
    """per K: list of modes (largest first) as (representative run, [runs], avg sim within)"""
    modes, k = {}, None
    rep_re = re.compile(r"^\s+(\S+) represents (\d+) runs?: runs? (.*?)(?:\. Avg sim within = ([0-9.eE+-]+))?\s*$")
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^_+K=(\d+)", line)
            if m:
                k = int(m.group(1))
                modes[k] = []
                continue
            m = rep_re.match(line)
            if m and k is not None:
                runs = [r.strip() for r in m.group(3).split(",")]
                sim = float(m.group(4)) if m.group(4) else np.nan
                modes[k].append((m.group(1), runs, sim))
    for k in modes:
        modes[k].sort(key=lambda x: -len(x[1]))  # stable sort: pong's order among ties
    return modes


def read_q(runs_dir, run_id):
    k, r = re.fullmatch(r"K(\d+)r(\d+)", run_id).groups()
    return np.loadtxt(os.path.join(runs_dir, f"K{k}", f"rep{r}", f"geno.{k}.Q"), ndmin=2)


def match(ref, q):
    """column order of q so that its first ref.shape[1] columns best match ref's columns
    (overlap sum(min) over individuals, Hungarian); unmatched columns of q follow"""
    sim = np.array([[np.minimum(ref[:, i], q[:, j]).sum() for j in range(q.shape[1])] for i in range(ref.shape[1])])
    rows, cols = linear_sum_assignment(-sim)
    order = list(cols[np.argsort(rows)])
    return order + [j for j in range(q.shape[1]) if j not in order]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs-dir", required=True)
    p.add_argument("--k", required=True)
    p.add_argument("--reps", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--fam", required=True)
    p.add_argument("--sample-pop", required=True)
    p.add_argument("--pong-summary", required=True)
    p.add_argument("--stable-min-frac", type=float, default=1.0)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    ks = [int(k) for k in a.k.split(",")]

    rows = []
    for k in ks:
        for r in range(1, a.reps + 1):
            cv, ll = parse_log(os.path.join(a.runs_dir, f"K{k}", f"rep{r}", "admixture.log"))
            rows.append({"K": k, "rep": r, "seed": a.seed + r - 1, "cv_error": cv, "loglik": ll})
    cv = pd.DataFrame(rows)
    cv.to_csv(os.path.join(a.outdir, "cv_error.tsv"), sep="\t", index=False, float_format="%.6g")

    modes = parse_pong(a.pong_summary)
    missing = [k for k in ks if k not in modes]
    if missing:
        sys.exit(f"ERROR: K={missing} not in {a.pong_summary}")
    summ = cv.groupby("K").cv_error.agg(cv_mean="mean", cv_sd="std", cv_min="min").reset_index()
    summ["n_runs"] = a.reps
    summ["n_modes"] = [len(modes[k]) for k in summ.K]
    summ["major_mode_runs"] = [len(modes[k][0][1]) for k in summ.K]
    summ["major_mode_frac"] = summ.major_mode_runs / summ.n_runs
    summ["major_mode_sim"] = [modes[k][0][2] for k in summ.K]
    summ["major_mode_rep_run"] = [modes[k][0][0] for k in summ.K]
    # K = 1 has a single cluster, so every run is identical
    summ["stable"] = np.where((summ.K == 1) | (summ.major_mode_frac >= a.stable_min_frac - 1e-9), "yes", "no")
    summ["selected"] = np.where(summ.K == summ.loc[summ.cv_mean.idxmin(), "K"], "yes", "")
    summ.to_csv(os.path.join(a.outdir, "k_summary.tsv"), sep="\t", index=False, float_format="%.6g")

    kcv = int(summ.loc[summ.cv_mean.idxmin(), "K"])
    s = summ.set_index("K")
    sel = {"K_cv": kcv, "cv_mean": s.loc[kcv, "cv_mean"], "stable": s.loc[kcv, "stable"],
           "n_modes": s.loc[kcv, "n_modes"], "major_mode_runs": s.loc[kcv, "major_mode_runs"], "n_runs": a.reps,
           "stable_min_frac": a.stable_min_frac}
    km1 = kcv - 1
    sel["K_minus_1"] = km1 if km1 in s.index else ""
    for col, src in (("K_minus_1_stable", "stable"), ("K_minus_1_n_modes", "n_modes"),
                     ("K_minus_1_major_mode_runs", "major_mode_runs")):
        sel[col] = s.loc[km1, src] if km1 in s.index else ""
    stable_le = [k for k in s.index if k <= kcv and s.loc[k, "stable"] == "yes"]
    sel["largest_stable_K_le_K_cv"] = max(stable_le) if stable_le else ""
    pd.DataFrame([sel]).to_csv(os.path.join(a.outdir, "selected_k.tsv"), sep="\t", index=False, float_format="%.6g")

    fam = pd.read_csv(a.fam, sep=r"\s+", header=None, dtype=str)
    pop_of = dict(pd.read_csv(a.sample_pop, sep="\t", dtype=str)[["vcf_id", "popid"]].values)
    ids = pd.DataFrame({"sample_id": fam[1], "popid": [pop_of[i] for i in fam[1]]})

    # major-mode representative run at each K, clusters aligned to the previous K
    aligned, prev = {}, None
    for k in ks:
        q = read_q(a.runs_dir, modes[k][0][0])
        if q.shape != (len(ids), k):
            sys.exit(f"ERROR: Q matrix for K={k} has shape {q.shape}, expected ({len(ids)}, {k})")
        if prev is not None:
            q = q[:, match(prev, q)]
        aligned[k] = prev = q
        out = pd.concat([ids, pd.DataFrame(q, columns=[f"cluster{j + 1}" for j in range(k)])], axis=1)
        out.to_csv(os.path.join(a.outdir, f"ancestry_K{k}.tsv"), sep="\t", index=False, float_format="%.6g")

    mode_rows = []
    for i, (rep, runs, sim) in enumerate(modes[kcv], start=1):
        q = read_q(a.runs_dir, rep)
        q = q[:, match(aligned[kcv], q)]
        d = pd.concat([ids, pd.DataFrame(q, columns=[f"cluster{j + 1}" for j in range(kcv)])], axis=1)
        d.insert(0, "avg_sim_within", sim)
        d.insert(0, "runs", ",".join(runs))
        d.insert(0, "n_runs", len(runs))
        d.insert(0, "rep_run", rep)
        d.insert(0, "mode", i)
        mode_rows.append(d)
    pd.concat(mode_rows).to_csv(os.path.join(a.outdir, "selected_k_modes.tsv"), sep="\t", index=False,
                                float_format="%.6g")

    print(summ.to_string(index=False))
    print(pd.DataFrame([sel]).T.to_string(header=False))
    if sel["stable"] != "yes":
        print(f"WARNING: K = {kcv} (lowest CV error) is not stable: pong found {sel['n_modes']} modes among "
              f"{a.reps} runs", file=sys.stderr)


if __name__ == "__main__":
    main()
