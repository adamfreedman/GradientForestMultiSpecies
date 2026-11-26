"""validate_inputs.py, vcf_to_popfreq.py and genotype_pca.py."""
import os

import numpy as np
import pandas as pd
import pytest

from conftest import TEST_DATA, run_script

SPA = os.path.join(TEST_DATA, "spA")


def validate(tmp_path, env_tif, popmap=None, latlongs=None, id_regex="^TEST_", snap=10, vcf=None):
    p = run_script("validate_inputs.py", "--vcf", vcf or os.path.join(SPA, "spA.vcf"),
                   "--popmap", popmap or os.path.join(SPA, "spA_id_to_pop.txt"),
                   "--latlongs", latlongs or os.path.join(SPA, "spA_pop_latlongs.csv"),
                   "--env-tif", env_tif, "--id-regex", id_regex, "--max-snap-km", snap,
                   "--out-report", tmp_path / "report.txt", "--out-sample-pop", tmp_path / "sample_pop.tsv",
                   "--out-unmapped", tmp_path / "unmapped.tsv", check=False)
    return p


def test_validate_test_species(tmp_path, env_tif):
    p = validate(tmp_path, env_tif)
    assert p.returncode == 0, p.stderr
    sp = pd.read_csv(tmp_path / "sample_pop.tsv", sep="\t", dtype=str)
    assert len(sp) == 12 * 10 + 3
    assert sp.vcf_id.str.startswith("TEST_").all()  # original VCF IDs are kept
    un = pd.read_csv(tmp_path / "unmapped.tsv", sep="\t", dtype=str)
    assert un.sample_id.tolist() == ["TEST_spA_extra_01"]
    assert un.reason.tolist() == ["in VCF, not in popmap"]
    report = (tmp_path / "report.txt").read_text()
    assert "will snap" in report and "# errors (0)" in report


def test_validate_without_id_regex_fails(tmp_path, env_tif):
    # without the regex no VCF sample matches the popmap, so every population is empty
    p = validate(tmp_path, env_tif, id_regex="")
    assert p.returncode != 0
    assert "has no individuals in the VCF" in p.stderr
    assert not (tmp_path / "sample_pop.tsv").exists()


def test_validate_site_too_far_from_data(tmp_path, env_tif):
    p = validate(tmp_path, env_tif, snap=2)
    assert p.returncode != 0 and "on NA cell" in p.stderr


def test_validate_site_outside_grid_and_duplicates(tmp_path, env_tif):
    ll = pd.read_csv(os.path.join(SPA, "spA_pop_latlongs.csv"), dtype=str)
    ll.loc[1, "long"] = "20.0"
    f = tmp_path / "ll.csv"
    ll.to_csv(f, index=False)
    p = validate(tmp_path, env_tif, latlongs=f)
    assert p.returncode != 0 and "outside the environmental grid" in p.stderr

    pm = pd.read_csv(os.path.join(SPA, "spA_id_to_pop.txt"), sep=r"\s+", dtype=str)
    pm = pd.concat([pm, pm.iloc[[0]]])
    f = tmp_path / "pm.txt"
    pm.to_csv(f, sep=" ", index=False)
    p = validate(tmp_path, env_tif, popmap=f)
    assert p.returncode != 0 and "duplicate sampleids" in p.stderr


def test_validate_population_without_coordinates(tmp_path, env_tif):
    # drop population 13 from latlongs: its 3 individuals become unmapped (warning, not error)
    ll = pd.read_csv(os.path.join(SPA, "spA_pop_latlongs.csv"), dtype=str)
    f = tmp_path / "ll.csv"
    ll[ll.popid != "13"].to_csv(f, index=False)
    p = validate(tmp_path, env_tif, latlongs=f)
    assert p.returncode == 0, p.stderr
    un = pd.read_csv(tmp_path / "unmapped.tsv", sep="\t", dtype=str)
    assert (un.reason == "population has no coordinates in latlongs").sum() == 3


VCF_HEAD = """##fileformat=VCFv4.2
##contig=<ID=c1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta1\ta2\tb1\tb2\tc1\tc2\td1
"""


def popfreq(tmp_path, body, min_ind=2, check=True):
    """run vcf_to_popfreq on a hand-written VCF with populations
    1 (a1, a2), 2 (b1, b2), 3 (c1, c2) and 4 (d1 only)"""
    vcf = tmp_path / "t.vcf"
    vcf.write_text(VCF_HEAD + "".join(f"c1\t{i + 1}\t{line}\n" for i, line in enumerate(body)))
    sp = tmp_path / "sp.tsv"
    sp.write_text("vcf_id\tpopid\n" + "".join(f"{s}\t{p}\n" for s, p in
                                              [("a1", 1), ("a2", 1), ("b1", 2), ("b2", 2), ("c1", 3),
                                               ("c2", 3), ("d1", 4)]))
    sites = tmp_path / "sites.tsv"
    sites.write_text("popid\tsitename\n1\tA\n2\tB\n3\tC\n4\tD\n")
    out = {k: tmp_path / f"{k}.tsv" for k in ("freq", "summary", "pops", "info")}
    p = run_script("vcf_to_popfreq.py", "--vcf", vcf, "--sample-pop", sp, "--env-sites", sites,
                   "--min-ind-per-pop", min_ind, "--out", out["freq"], "--out-summary", out["summary"],
                   "--out-pops", out["pops"], "--out-snp-info", out["info"], check=check)
    return p, out


def test_popfreq_frequencies_and_filters(tmp_path):
    gt = "\t".join
    body = [
        "s1\tA\tG\t.\tPASS\t.\tGT\t" + gt(["0/0", "0/0", "0/1", "0/1", "1/1", "1/1", "1/1"]),  # cline
        "s2\tA\tG\t.\tPASS\t.\tGT\t" + gt(["1/1", "1/1", "1/1", "0/1", "0/0", "1/1", "0/0"]),  # ALT is major
        ".\tA\tG,T\t.\tPASS\t.\tGT\t" + gt(["0/1"] * 7),                                       # multiallelic
        "s4\tA\tG\t.\tPASS\t.\tGT\t" + gt(["./.", "./.", "0/1", "0/1", "0/0", "0/0", "0/0"]),  # pop 1 uncalled
        "s5\tA\tG\t.\tPASS\t.\tGT\t" + gt(["0/1"] * 7),                                        # no variation
        "s6\tA\tG\t.\tPASS\t.\tGT\t" + gt(["0/1", "./.", "0/0", "0/0", "0/0", "1/1", "0/0"]),  # missing call
        "s1\tA\tG\t.\tPASS\t.\tGT\t" + gt(["0/0", "0/1", "0/0", "0/0", "0/0", "0/0", "0/0"]),  # duplicate ID
        ".\tA\tC\t.\tPASS\t.\tGT\t" + gt(["0/0", "0/0", "0/0", "0/1", "0/0", "0/0", "0/0"]),   # no ID
    ]
    p, out = popfreq(tmp_path, body)
    assert "dropping populations with < 2 individuals: ['4']" in p.stderr
    f = pd.read_csv(out["freq"], sep="\t", dtype={"popid": str}).set_index("popid")
    assert f.index.tolist() == ["1", "2", "3"]
    assert f.columns.tolist() == ["snp_s1", "snp_s2", "snp_s6", "snp_s1_2", "c1_8"]
    assert f["snp_s1"].tolist() == [0.0, 0.5, 1.0]  # global frequency 0.5: not flipped
    assert f["snp_s2"].tolist() == [0.0, 0.25, 0.5]  # polarized to the minor allele (A)
    assert f["snp_s6"].tolist() == [0.5, 0.0, 0.5]   # frequency over called genotypes only
    summ = dict(pd.read_csv(out["summary"], sep="\t").values)
    assert summ == {"total": 8, "not_biallelic_snp": 1, "no_calls_in_a_pop": 1,
                    "no_among_pop_variation": 1, "kept": 5}
    pops = pd.read_csv(out["pops"], sep="\t", dtype={"popid": str})
    assert pops.n_ind.tolist() == [2, 2, 2]
    info = pd.read_csv(out["info"], sep="\t", dtype=str)
    assert info.snp.tolist() == f.columns.tolist()
    assert info.VCF_ID.tolist()[-1] == "."


def test_popfreq_too_few_populations(tmp_path):
    p, _ = popfreq(tmp_path, ["s1\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(["0/1"] * 7)], min_ind=3, check=False)
    assert p.returncode != 0 and "fewer than 3 populations" in p.stderr


def test_popfreq_no_snps_left(tmp_path):
    p, _ = popfreq(tmp_path, ["s1\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(["0/1"] * 7)], check=False)
    assert p.returncode != 0 and "no SNPs left" in p.stderr


@pytest.mark.parametrize("vcf", ["spA.vcf", "spB.vcf.gz"])
def test_popfreq_on_test_species(tmp_path, env_tif, vcf):
    sp_name = vcf.split(".")[0]
    d = os.path.join(TEST_DATA, sp_name)
    regex = "^TEST_" if sp_name == "spA" else ""
    p = validate(tmp_path, env_tif, popmap=os.path.join(d, f"{sp_name}_id_to_pop.txt"),
                 latlongs=os.path.join(d, f"{sp_name}_pop_latlongs.csv"), id_regex=regex,
                 vcf=os.path.join(d, vcf))
    assert p.returncode == 0, p.stderr
    run_script("extract_env.py", "--latlongs", os.path.join(d, f"{sp_name}_pop_latlongs.csv"),
               "--sample-pop", tmp_path / "sample_pop.tsv", "--env-tif", env_tif, "--max-snap-km", 10,
               "--out", tmp_path / "env_sites.tsv")
    out = {k: tmp_path / f"{k}.tsv" for k in ("freq", "summary", "pops", "info")}
    run_script("vcf_to_popfreq.py", "--vcf", os.path.join(d, vcf), "--sample-pop", tmp_path / "sample_pop.tsv",
               "--env-sites", tmp_path / "env_sites.tsv", "--min-ind-per-pop", 4, "--out", out["freq"],
               "--out-summary", out["summary"], "--out-pops", out["pops"], "--out-snp-info", out["info"])
    f = pd.read_csv(out["freq"], sep="\t", dtype={"popid": str})
    summ = dict(pd.read_csv(out["summary"], sep="\t").values)
    if sp_name == "spA":
        assert len(f) == 12 and "13" not in f.popid.tolist()  # 3 individuals < 4
        assert summ["total"] == 402 and summ["not_biallelic_snp"] == 1
        assert summ["no_among_pop_variation"] >= 1
    else:
        assert len(f) == 10 and summ["total"] == 250
    assert f.drop(columns="popid").stack().between(0, 1).all()

    # genotype PCA on the same individuals
    pca = {k: tmp_path / f"pca_{k}.tsv" for k in ("variance", "scores", "info")}
    run_script("genotype_pca.py", "--vcf", os.path.join(d, vcf), "--sample-pop", tmp_path / "sample_pop.tsv",
               "--env-sites", tmp_path / "env_sites.tsv", "--pop-counts", out["pops"],
               "--out-variance", pca["variance"], "--out-scores", pca["scores"], "--out-info", pca["info"])
    sc = pd.read_csv(pca["scores"], sep="\t", dtype={"popid": str, "sample_id": str})
    va = pd.read_csv(pca["variance"], sep="\t")
    n = 123 if sp_name == "spA" else 80
    assert len(sc) == n and sc.sample_id.is_unique  # every mapped individual, incl. small populations
    assert set(sc.popid[sc.in_gf_model == "no"]) == ({"13"} if sp_name == "spA" else set())
    assert len(va) <= n - 1 and va.pct_variance.is_monotonic_decreasing
    assert va.pct_variance.sum() == pytest.approx(100, abs=0.01)
    assert sc.columns.tolist()[5:] == va.PC.tolist()
    # the simulated populations differ in allele frequency: PC1 separates them
    by_pop = sc.groupby("popid").PC1
    assert by_pop.mean().var() > 5 * by_pop.var().mean()


def test_genotype_pca_matches_direct_svd(tmp_path):
    """scores and variance shares equal an SVD of the standardized, mean-imputed matrix"""
    rng = np.random.default_rng(3)
    samples = ["a1", "a2", "b1", "b2", "c1", "c2", "d1"]
    G = rng.integers(0, 3, size=(40, 7)).astype(float)
    G[rng.random(G.shape) < 0.1] = np.nan
    G[5] = 1   # every individual heterozygous: no variation, skipped
    G[6] = np.nan  # no calls: skipped
    code = {0.0: "0/0", 1.0: "0/1", 2.0: "1/1"}
    body = ["\t".join(["c1", str(i + 1), f"s{i}", "A", "G", ".", "PASS", ".", "GT"] +
                      [code.get(g, "./.") for g in row]) for i, row in enumerate(G)]
    body.append("\t".join(["c1", "99", "m", "A", "G,T", ".", "PASS", ".", "GT"] + ["0/1"] * 7))  # multiallelic
    vcf = tmp_path / "t.vcf"
    vcf.write_text(VCF_HEAD + "\n".join(body) + "\n")
    sp = tmp_path / "sp.tsv"
    sp.write_text("vcf_id\tpopid\n" + "".join(f"{s}\t{s[0]}\n" for s in samples))
    sites = tmp_path / "sites.tsv"
    sites.write_text("popid\tsitename\na\tA\nb\tB\nc\tC\nd\tD\n")
    pops = tmp_path / "pops.tsv"
    pops.write_text("popid\tn_ind\na\t2\nb\t2\nc\t2\n")
    out = {k: tmp_path / f"{k}.tsv" for k in ("variance", "scores", "info")}
    run_script("genotype_pca.py", "--vcf", vcf, "--sample-pop", sp, "--env-sites", sites, "--pop-counts", pops,
               "--out-variance", out["variance"], "--out-scores", out["scores"], "--out-info", out["info"])

    keep = [i for i in range(len(G)) if i not in (5, 6)]
    X = G[keep].T  # individuals x SNPs
    f = np.nanmean(X, axis=0) / 2
    Z = np.nan_to_num((X - 2 * f) / np.sqrt(2 * f * (1 - f)))
    U, S, _ = np.linalg.svd(Z / np.sqrt(Z.shape[1]), full_matrices=False)
    expect = U * S
    pct = 100 * S ** 2 / (S ** 2).sum()

    sc = pd.read_csv(out["scores"], sep="\t", dtype=str)
    va = pd.read_csv(out["variance"], sep="\t")
    k = len(va)
    assert k == 6  # 7 individuals: at most 6 PCs after centring
    got = sc[[f"PC{i + 1}" for i in range(k)]].astype(float).to_numpy()
    assert np.allclose(np.abs(got), np.abs(expect[:, :k]), atol=1e-4)
    assert np.allclose(va.pct_variance, pct[:k], atol=1e-4)
    assert (got[np.abs(got).argmax(axis=0), range(k)] > 0).all()  # sign convention
    assert sc.in_gf_model.tolist() == ["yes"] * 6 + ["no"]
    assert sc.sitename.tolist() == ["A", "A", "B", "B", "C", "C", "D"]
    assert np.allclose(sc.missing_frac.astype(float), np.isnan(X).mean(axis=1))
    info = dict(pd.read_csv(out["info"], sep="\t").values)
    assert info == {"n_individuals": 7, "snps_total": 41, "snps_not_biallelic_snp": 1, "snps_no_calls": 1,
                    "snps_invariant": 1, "snps_used": 38}
