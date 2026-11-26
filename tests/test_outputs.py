"""Checks on the outputs of a full workflow run on the test data.

Run the workflow first (from the repository root):
    snakemake --use-conda --cores 4 --configfile .test/config/config.yaml
then:
    pytest tests/test_outputs.py
By default this checks .test/results. Set GF_TEST_OUTDIRS to a comma-separated
list of output directories to check alternative runs as well, e.g.
    GF_TEST_OUTDIRS=.test/results,.test/results_importance_floor0.07,.test/results_all_predictors
Directories that do not exist are skipped.
"""
import itertools
import os

import numpy as np
import pandas as pd
import pytest
import rasterio

from conftest import ROOT

OUTDIRS = [d.strip() for d in os.environ.get("GF_TEST_OUTDIRS", ".test/results").split(",") if d.strip()]
SPECIES = {"spA": 12, "spB": 10}  # populations left after min_ind_per_pop
N_RANDOM = 6
PLOTS = ["map_turnover_rgb", "biplot_turnover_pca", "importance_r2weighted", "importance_accuracy",
         "random_vs_real"]
# the variables that drive the simulated allele frequencies (see .test/make_test_data.py)
DRIVERS = {"spA": {"bio1", "bio2"}, "spB": {"bio3", "bio4"}}


def existing(outdir):
    path = os.path.join(ROOT, outdir)
    if not os.path.isdir(path):
        pytest.skip(f"{outdir} not found; run the workflow first")
    return path


def kind(outdir):
    if outdir.endswith("_all_predictors"):
        return "all"
    if "_importance_floor" in outdir:
        return "floor"
    return "default"


def read(path, **kw):
    return pd.read_csv(path, sep="\t", **kw)


SP_CASES = [(o, s) for o in OUTDIRS for s in SPECIES]
MODEL_CASES = [(o, s, m) for o, s in SP_CASES for m in ("", "refit_sigvars")]


@pytest.mark.parametrize("outdir,sp", SP_CASES)
def test_inputs_and_genotypes(outdir, sp):
    d = os.path.join(existing(outdir), sp)
    freq = read(f"{d}/genotypes/popfreq.tsv", dtype={"popid": str})
    assert len(freq) == SPECIES[sp]
    assert freq.drop(columns="popid").stack().between(0, 1).all()
    info = read(f"{d}/genotypes/snp_info.tsv")
    assert info.snp.tolist() == freq.columns[1:].tolist()
    sites = read(f"{d}/env_sites.tsv", dtype={"popid": str})
    un = read(f"{d}/qc/unmapped_samples.tsv")
    if sp == "spA":
        assert "13" not in freq.popid.tolist()  # 3 individuals < min_ind_per_pop 4
        assert un.sample_id.tolist() == ["TEST_spA_extra_01"]
        assert sites.set_index("popid").loc["1", "snapped_km"] > 0
    else:
        assert un.empty
    assert "# errors (0)" in open(f"{d}/qc/validation.txt").read()


def selected_vars(d):
    return [v for v in open(f"{d}/predictors/selected.txt").read().split() if v]


@pytest.mark.parametrize("outdir,sp", SP_CASES)
def test_predictor_selection(outdir, sp):
    d = os.path.join(existing(outdir), sp)
    rep = read(f"{d}/predictors/report.tsv").set_index("variable")
    sel = selected_vars(d)
    assert sorted(sel) == sorted(rep.index[rep.status == "kept"])
    assert rep.loc["bio5", "status"] == "dropped"
    assert rep.loc["bio5", "reason"] == "constant across sites"
    sites = read(f"{d}/env_sites.tsv")
    if kind(outdir) == "all":
        assert sorted(sel) == ["bio1", "bio2", "bio3", "bio4", "bio6"]
        return
    # selection: no two kept variables correlated above 0.8 at the sites
    r = sites[sel].corr().abs()
    for a, b in itertools.combinations(sel, 2):
        assert r.loc[a, b] <= 0.8, (a, b)
    dropped = rep[(rep.status == "dropped") & (rep.index != "bio5")]
    floor = kind(outdir) == "floor"
    by_floor = dropped.reason.str.contains("importance floor")
    assert all(dropped.reason.str.startswith("|r|=") | (floor & by_floor))
    assert "screen_share" in rep.columns
    if floor:
        # the floor can remove both members of the correlated pair bio1/bio2
        assert len({"bio1", "bio2"} & set(sel)) <= 1
        f = float(outdir.split("_importance_floor")[1])
        assert (rep.loc[sel, "screen_share"] >= f).all()
        assert (rep.loc[by_floor[by_floor].index, "screen_share"] < f).all()
    else:
        assert len({"bio1", "bio2"} & set(sel)) == 1


def test_floor_drops_variables():
    """the test floor (0.07) is set so that it removes at least one variable"""
    floor_dirs = [o for o in OUTDIRS if kind(o) == "floor"]
    if not floor_dirs:
        pytest.skip("no importance-floor run requested in GF_TEST_OUTDIRS")
    n = 0
    for sp in SPECIES:
        rep = read(os.path.join(existing(floor_dirs[0]), sp, "predictors/report.tsv"))
        n += rep.reason.str.contains("importance floor", na=False).sum()
    assert n > 0


@pytest.mark.parametrize("outdir,sp", SP_CASES)
def test_stability(outdir, sp):
    d = os.path.join(existing(outdir), sp)
    if kind(outdir) == "all":
        assert not os.path.exists(f"{d}/predictors/stability_variables.tsv")
        return
    v = read(f"{d}/predictors/stability_variables.tsv")
    s = read(f"{d}/predictors/stability_sites.tsv")
    assert len(s) == SPECIES[sp]  # one row per population with genotype data
    assert set(s.same_as_full) <= {"yes", "no"}
    assert v.kept_in.between(0, SPECIES[sp]).all()
    assert set(selected_vars(d)) <= set(v.variable)


@pytest.mark.parametrize("outdir,sp,sub", MODEL_CASES)
def test_model_outputs(outdir, sp, sub):
    d = os.path.join(existing(outdir), sp, sub)
    preds = (open(f"{d}/predictors.txt").read().split() if sub
             else selected_vars(os.path.join(existing(outdir), sp)))
    if sub:  # alpha 1 in the test config: every predictor passes
        assert sorted(preds) == sorted(selected_vars(os.path.join(existing(outdir), sp)))
    for p in PLOTS:
        for ext in ("pdf", "png"):
            assert os.path.getsize(f"{d}/plots/{p}.{ext}") > 1000
    assert os.path.getsize(f"{d}/plots/gf_diagnostics.pdf") > 1000

    summ = read(f"{d}/gf/summary.tsv").iloc[0]
    assert summ.n_sites == SPECIES[sp] and summ.n_predictors == len(preds)
    assert summ.n_pos_rsq > 0 and 0 < summ.mean_rsq <= 1
    assert summ.ntree in (50, 100)
    assert read(f"{d}/gf/ntree_tuning.tsv").shape[0] >= 1

    imp = read(f"{d}/gf/importance.tsv")
    assert sorted(imp.variable) == sorted(preds)

    null = read(f"{d}/gf/random_null.tsv")
    assert len(null) == N_RANDOM and sorted(null.rep) == list(range(1, N_RANDOM + 1))
    pv = read(f"{d}/gf/random_pvalues.tsv").set_index("statistic")
    for stat in ["n_pos_rsq", "mean_rsq", "sum_rsq"] + [f"imp_{v}" for v in preds] + [f"acc_{v}" for v in preds]:
        assert stat in pv.index, stat
    assert pv.p_empirical.between(1 / (N_RANDOM + 1) - 1e-9, 1).all()
    assert (pv.n_random == N_RANDOM).all()

    snp = read(f"{d}/gf/snp_pvalues.tsv")
    n_snps = read(os.path.join(existing(outdir), sp, "genotypes/popfreq.tsv"), nrows=0).shape[1] - 1
    assert len(snp) == n_snps
    tested = snp.rsq > 0
    assert tested.sum() == summ.n_pos_rsq
    assert snp.loc[tested, "p_bh"].between(0, 1).all()
    assert snp.loc[~tested, ["p_empirical", "p_bh"]].isna().all().all()  # only SNPs with R^2 > 0 are tested

    with rasterio.open(f"{d}/gf/turnover_stack.tif") as src:
        assert src.count == len(preds)
    with rasterio.open(f"{d}/gf/turnover_rgb.tif") as src:
        assert src.count == 3
    assert len(read(f"{d}/gf/site_pcs.tsv")) == SPECIES[sp]
    assert sorted(read(f"{d}/gf/pca_loadings.tsv").variable) == sorted(preds)

    html = open(f"{d}/{'qc_refit' if sub else 'qc'}.html").read()
    assert f"<h1>{sp}" in html and "How to read this report" in html
    assert html.count("data:image/png;base64,") >= 5
    if kind(outdir) == "floor" and not sub:
        assert "Importance floor" in html


@pytest.mark.parametrize("outdir,sp", SP_CASES)
def test_signal_recovered(outdir, sp):
    """the variable that drives the simulated allele frequencies ranks first"""
    d = os.path.join(existing(outdir), sp)
    imp = read(f"{d}/gf/importance.tsv").sort_values("importance", ascending=False)
    assert imp.variable.iloc[0] in DRIVERS[sp]
    pv = read(f"{d}/gf/random_pvalues.tsv").set_index("statistic")
    assert pv.loc["n_pos_rsq", "observed"] > pv.loc["n_pos_rsq", "null_mean"]
    assert np.isfinite(pv.loc["mean_rsq", "observed"])
