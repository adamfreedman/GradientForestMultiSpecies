"""admixture_summary.py: parsing pong's summary, cluster alignment across K, and the
choice of K with its stability, on hand-made Q matrices and pong output."""
import os

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scipy")
from admixture_summary import match, parse_pong  # noqa: E402
from conftest import run_script  # noqa: E402

PONG_SUMMARY = """
Results summary
===============

____________________________________________K=2
Number of runs: 3
Major mode: K2r2

Representative runs: K2r2
\tK2r2 represents 3 runs: runs K2r2, K2r3, K2r1. Avg sim within = 0.999999

____________________________________________K=3
Number of runs: 3
Major mode: K3r2

Representative runs: K3r1, K3r2
\tK3r1 represents 1 run: run K3r1
\tK3r2 represents 2 runs: runs K3r2, K3r3. Avg sim within = 0.990000
Avg sim between modes = 0.811587
"""


def test_parse_pong(tmp_path):
    f = tmp_path / "result_summary.txt"
    f.write_text(PONG_SUMMARY)
    m = parse_pong(f)
    assert m[2] == [("K2r2", ["K2r2", "K2r3", "K2r1"], pytest.approx(0.999999))]
    # largest mode first, whatever pong's order
    assert [x[0] for x in m[3]] == ["K3r2", "K3r1"]
    assert m[3][0][1] == ["K3r2", "K3r3"] and np.isnan(m[3][1][2])


def test_match_aligns_columns():
    ref = np.array([[1.0, 0.0], [0.9, 0.1], [0.0, 1.0], [0.1, 0.9]])
    # same two clusters in the other order, plus a new third cluster
    q = np.column_stack([ref[:, 1] * 0.8, np.full(4, 0.2), ref[:, 0] * 0.8])
    assert match(ref, q) == [2, 0, 1]
    assert match(ref, ref[:, ::-1]) == [1, 0]


def write_run(d, k, r, q, cv):
    run = d / f"K{k}" / f"rep{r}"
    run.mkdir(parents=True)
    np.savetxt(run / f"geno.{k}.Q", q, fmt="%.6f")
    (run / "admixture.log").write_text(f"Loglikelihood: -{1000 - k}.5\nCV error (K={k}): {cv}\n")


def test_summary_selects_k_and_flags_instability(tmp_path):
    rng = np.random.default_rng(1)
    n = 12
    pop = np.repeat(["1", "2", "10"], 4)
    runs = tmp_path / "runs"
    base3 = np.zeros((n, 3))
    base3[np.arange(n), np.repeat([0, 1, 2], 4)] = 1
    base3 = 0.9 * base3 + 0.1 / 3
    for r in (1, 2, 3):
        write_run(runs, 2, r, np.column_stack([base3[:, 0] + base3[:, 1], base3[:, 2]]), 0.50 + r / 1000)
        # K = 3 has the lowest CV error; rep 1 has its clusters in a different order
        q3 = base3[:, [2, 0, 1]] if r == 1 else base3
        write_run(runs, 3, r, q3 + rng.uniform(0, 1e-3, q3.shape), 0.40 + r / 1000)
    fam = tmp_path / "geno.fam"
    pd.DataFrame({"f": 0, "i": [f"s{i}" for i in range(n)], "a": 0, "b": 0, "c": 0, "d": -9}).to_csv(
        fam, sep=" ", header=False, index=False)
    sp = tmp_path / "sample_pop.tsv"
    pd.DataFrame({"vcf_id": [f"s{i}" for i in range(n)], "popid": pop}).to_csv(sp, sep="\t", index=False)
    pong = tmp_path / "result_summary.txt"
    pong.write_text(PONG_SUMMARY)
    out = tmp_path / "out"
    out.mkdir()
    run_script("admixture_summary.py", "--runs-dir", runs, "--k", "2,3", "--reps", 3, "--seed", 5,
               "--fam", fam, "--sample-pop", sp, "--pong-summary", pong, "--stable-min-frac", 1.0, "--outdir", out)

    cv = pd.read_csv(out / "cv_error.tsv", sep="\t")
    assert len(cv) == 6 and cv.seed.tolist() == [5, 6, 7] * 2
    ks = pd.read_csv(out / "k_summary.tsv", sep="\t").set_index("K")
    assert ks.loc[2, "stable"] == "yes" and ks.loc[3, "stable"] == "no"
    assert ks.loc[3, "n_modes"] == 2 and ks.loc[3, "major_mode_runs"] == 2
    sel = pd.read_csv(out / "selected_k.tsv", sep="\t").iloc[0]
    assert sel.K_cv == 3 and sel.stable == "no" and sel.K_minus_1 == 2 and sel.K_minus_1_stable == "yes"
    assert sel.largest_stable_K_le_K_cv == 2

    a2 = pd.read_csv(out / "ancestry_K2.tsv", sep="\t", dtype={"popid": str})
    a3 = pd.read_csv(out / "ancestry_K3.tsv", sep="\t", dtype={"popid": str})
    assert a3.sample_id.tolist() == [f"s{i}" for i in range(n)] and a3.popid.tolist() == list(pop)
    assert np.allclose(a3.filter(like="cluster").sum(axis=1), 1, atol=1e-2)
    # cluster 1 at K = 2 (populations 1 + 2) splits at K = 3; the population-10 cluster keeps
    # its column number (2), and the new cluster is the last one
    assert a2.cluster2[pop == "10"].min() > 0.8
    assert a3.cluster2[pop == "10"].min() > 0.8
    assert a3.cluster3[pop == "2"].min() > 0.8 or a3.cluster3[pop == "1"].min() > 0.8
    modes = pd.read_csv(out / "selected_k_modes.tsv", sep="\t")
    assert modes["mode"].unique().tolist() == [1, 2] and modes.groupby("mode").n_runs.first().tolist() == [2, 1]
    # the odd run, once aligned, matches the major mode
    m1 = modes[modes["mode"] == 1].filter(like="cluster").to_numpy()
    m2 = modes[modes["mode"] == 2].filter(like="cluster").to_numpy()
    assert np.abs(m1 - m2).max() < 0.01
