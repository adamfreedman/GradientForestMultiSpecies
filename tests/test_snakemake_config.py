"""Snakemake dry runs with the test configuration: the DAG builds, the
alternative runs write to their own directories, and bad settings are rejected
before anything runs."""
import shutil
import subprocess

import pytest

from conftest import ROOT

pytestmark = pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed")
TEST_CONFIG = ".test/config/config.yaml"


def dry_run(*config):
    args = ["snakemake", "-n", "--quiet", "rules", "--configfile", TEST_CONFIG]
    if config:
        args += ["--config", *config]
    return subprocess.run(args, capture_output=True, text=True, cwd=ROOT)


def dag_files(*config):
    """every output path in the dry-run plan (from --summary)"""
    args = ["snakemake", "--summary", "--configfile", TEST_CONFIG]
    if config:
        args += ["--config", *config]
    p = subprocess.run(args, capture_output=True, text=True, cwd=ROOT)
    assert p.returncode == 0, p.stderr
    return [line.split("\t")[0] for line in p.stdout.splitlines()[1:] if line.strip()]


def test_default_dag():
    p = dry_run()
    assert p.returncode == 0, p.stderr
    files = dag_files()
    for f in (".test/results/spA/qc.html", ".test/results/spB/refit_sigvars/qc_refit.html",
              ".test/results/spA/predictors/stability_variables.tsv", ".test/results/env/env_stack.tif"):
        assert f in files
    # one leave-one-out screening fit per population in each latlongs file
    assert sum("spA/predictors/leave_one_out/site" in f for f in files) == 13
    assert sum("spB/predictors/leave_one_out/site" in f for f in files) == 10


@pytest.mark.parametrize("config,outdir", [
    (["use_all_predictors=true"], ".test/results_all_predictors"),
    (["importance_floor=0.07"], ".test/results_importance_floor0.07"),
])
def test_alternative_runs_use_own_directory(config, outdir):
    files = dag_files(*config)
    model_files = [f for f in files if "/env/" not in f and not f.endswith(".gf_installed")
                   and "/logs/" not in f]
    assert model_files and all(f.startswith(outdir + "/") for f in model_files)
    # shared one-time setup stays under the base outdir
    assert ".test/results/env/env_stack.tif" in files
    if "use_all_predictors=true" in config:
        assert not any("leave_one_out" in f for f in files)  # no selection, no stability check


@pytest.mark.parametrize("config,message", [
    (["importance_floor=1.5"], "importance_floor must be a fraction"),
    (["importance_floor=0.01", "use_all_predictors=true"], "cannot be combined"),
    (["predictors={'rank_by': 'bogus'}"], "rank_by must be"),
    (["gf={'refit_significant': {'subdir': 'a/b'}}"], "plain directory name"),
])
def test_bad_settings_are_rejected(config, message):
    p = dry_run(*config)
    assert p.returncode != 0
    assert message in p.stdout + p.stderr
