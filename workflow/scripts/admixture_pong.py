#!/usr/bin/env python
"""Run pong on the ADMIXTURE runs (no web server) to group the replicate runs at each K
into modes.

Writes pong's inputs next to the filemap (filemap.txt, ind2pop.txt, pop_order.txt) and
pong's output to --outdir (result_summary.txt, params_used.txt). pong's verbose mode
is not used: it writes one file per pair of runs, ~5,000 files per K at 100 runs.
To browse the runs in pong's interactive viewer, run from the directory of the filemap:
    pong -m filemap.txt -i ind2pop.txt -n pop_order.txt -o <new dir>
Run IDs are K<k>r<rep>.
"""
import argparse
import os
import subprocess
import sys

import pandas as pd


def pop_order(popids):
    """population order used by the PCA plots: numeric IDs as numbers, otherwise sorted"""
    u = list(dict.fromkeys(popids))
    try:
        return sorted(u, key=float)
    except ValueError:
        return sorted(u)


def pong_label(x):
    return "_".join(str(x).split())  # pong labels cannot contain whitespace


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs-dir", required=True)
    p.add_argument("--k", required=True, help="comma-separated K values")
    p.add_argument("--reps", type=int, required=True)
    p.add_argument("--fam", required=True)
    p.add_argument("--sample-pop", required=True)
    p.add_argument("--sim-threshold", type=float, default=0.97)
    p.add_argument("--filemap", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()

    ks = [int(k) for k in a.k.split(",")]
    base = os.path.dirname(os.path.abspath(a.filemap))
    os.makedirs(base, exist_ok=True)
    with open(a.filemap, "w") as fh:
        for k in ks:
            for r in range(1, a.reps + 1):
                q = os.path.join(a.runs_dir, f"K{k}", f"rep{r}", f"geno.{k}.Q")
                if not os.path.exists(q):
                    sys.exit(f"ERROR: missing {q}")
                fh.write(f"K{k}r{r}\t{k}\t{os.path.relpath(os.path.abspath(q), base)}\n")

    fam = pd.read_csv(a.fam, sep=r"\s+", header=None, dtype=str)
    pop_of = dict(pd.read_csv(a.sample_pop, sep="\t", dtype=str)[["vcf_id", "popid"]].values)
    pops = [pong_label(pop_of[i]) for i in fam[1]]  # Q-matrix rows are in .fam order
    ind2pop = os.path.join(base, "ind2pop.txt")
    order = os.path.join(base, "pop_order.txt")
    pd.Series(pops).to_csv(ind2pop, index=False, header=False)
    pd.Series(pop_order(pops)).to_csv(order, index=False, header=False)

    cmd = ["pong", "-m", a.filemap, "-i", ind2pop, "-n", order, "-o", a.outdir, "-s", str(a.sim_threshold),
           "-f", "-g", "--disable_server"]
    print(" ".join(cmd), flush=True)
    # -f overwrites, -g answers pong's only other prompt (greedy matching when the
    # modes are not disjoint); stdin is closed so a new prompt fails instead of hanging
    res = subprocess.run(cmd, stdin=subprocess.DEVNULL)
    if res.returncode != 0 or not os.path.exists(os.path.join(a.outdir, "result_summary.txt")):
        sys.exit(f"ERROR: pong failed (exit code {res.returncode})")


if __name__ == "__main__":
    main()
