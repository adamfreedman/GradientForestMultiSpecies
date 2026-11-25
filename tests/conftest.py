"""Shared fixtures for the Python tests. Run from the repository root: pytest tests"""
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "workflow", "scripts")
TEST_DATA = os.path.join(ROOT, ".test", "data")
sys.path.insert(0, SCRIPTS)


def run_script(name, *args, check=True):
    """run workflow/scripts/<name> with the current Python; returns CompletedProcess"""
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, name), *map(str, args)],
                       capture_output=True, text=True, cwd=ROOT)
    if check and p.returncode != 0:
        raise AssertionError(f"{name} failed ({p.returncode}):\n{p.stdout}\n{p.stderr}")
    return p


@pytest.fixture(scope="session")
def env_tif(tmp_path_factory):
    """the test environmental stack, built from .test/data like the workflow does"""
    d = tmp_path_factory.mktemp("env")
    run_script("bbox_to_vector.py", "--vector", os.path.join(TEST_DATA, "study_area.geojson"),
               "--geojson", d / "area.geojson", "--bbox", d / "bbox.tsv")
    run_script("prep_env_stack.py", "--env-dir", os.path.join(TEST_DATA, "env"), "--pattern", "*.asc",
               "--name-regex", r"bio_?(\d+)", "--name-template", "bio{0}", "--clip", d / "area.geojson",
               "--aggregate", 1, "--out-tif", d / "env_stack.tif", "--out-layers", d / "layers.txt")
    return d / "env_stack.tif"
