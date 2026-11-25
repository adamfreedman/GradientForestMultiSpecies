"""select_sig_predictors.py: predictors for the refit."""
import pandas as pd

from conftest import run_script

PVALS = """statistic\tobserved\tp_empirical
n_pos_rsq\t100\t0.01
imp_bio1\t0.5\t0.01
imp_bio3\t0.2\t0.20
imp_bio4\t0.1\t0.60
acc_bio1\t0.01\t0.02
acc_bio3\t0.005\t0.04
acc_bio4\t0.001\t0.70
"""


def run(tmp_path, metrics, alpha=0.05, min_n=2):
    f = tmp_path / "p.tsv"
    f.write_text(PVALS)
    return run_script("select_sig_predictors.py", "--pvals", f, "--alpha", alpha, "--metrics", metrics,
                      "--min-predictors", min_n, "--out-list", tmp_path / "list.txt",
                      "--out-report", tmp_path / "rep.tsv", check=False)


def test_any_metric_passes(tmp_path):
    p = run(tmp_path, "importance,accuracy")
    assert p.returncode == 0, p.stderr
    assert (tmp_path / "list.txt").read_text().split() == ["bio1", "bio3"]  # bio3 via accuracy
    rep = pd.read_csv(tmp_path / "rep.tsv", sep="\t")
    assert rep.set_index("variable").loc["bio3", "reason"] == "accuracy p=0.04"
    assert "n_pos_rsq" not in rep.variable.tolist()


def test_too_few_predictors(tmp_path):
    p = run(tmp_path, "importance")
    assert p.returncode != 0 and "only 1 predictor(s) pass" in p.stderr
    assert not (tmp_path / "list.txt").exists()
    assert (tmp_path / "rep.tsv").exists()  # the report is still written


def test_unknown_metric(tmp_path):
    p = run(tmp_path, "importance,bogus")
    assert p.returncode != 0 and "unknown metric" in p.stderr
