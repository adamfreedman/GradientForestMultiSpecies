"""bbox_to_vector.py, prep_env_stack.py, extract_env.py and get_boundaries.py."""
import json
import os

import numpy as np
import pandas as pd
import pytest
import rasterio

from conftest import TEST_DATA, run_script
from prep_env_stack import block_mean, layer_name


def test_bbox_from_extent(tmp_path):
    run_script("bbox_to_vector.py", "--extent", 10, 0, 14, 3, "--geojson", tmp_path / "a.geojson",
               "--bbox", tmp_path / "b.tsv")
    bb = pd.read_csv(tmp_path / "b.tsv", sep="\t")
    assert bb.iloc[0].tolist() == [10, 0, 14, 3]
    gj = json.loads((tmp_path / "a.geojson").read_text())
    assert len(gj["features"]) == 1


def test_bbox_rejects_inverted_extent(tmp_path):
    p = run_script("bbox_to_vector.py", "--extent", 14, 0, 10, 3, "--geojson", tmp_path / "a.geojson",
                   check=False)
    assert p.returncode != 0 and "min < max" in p.stderr


@pytest.mark.parametrize("fname,expected", [
    ("PIREworldclim-bio_12.asc", ("bio12", 12)),
    ("PIREworldclim-bio1.asc", ("bio1", 1)),
    ("test-BIO_3.tif", ("bio3", 3)),
])
def test_layer_name(fname, expected):
    assert layer_name(fname, r"bio_?(\d+)", "bio{0}") == expected


def test_layer_name_no_match():
    with pytest.raises(SystemExit):
        layer_name("elevation.asc", r"bio_?(\d+)", "bio{0}")


def test_block_mean_ignores_nan():
    a = np.array([[1, 2, 3, 4], [3, np.nan, 5, 6], [np.nan, np.nan, 1, 1], [np.nan, np.nan, 1, 1]], float)
    m = block_mean(a, 2)
    assert m[0, 0] == pytest.approx(2.0)
    assert m[0, 1] == pytest.approx(4.5)
    assert np.isnan(m[1, 0])
    assert m[1, 1] == 1


def test_env_stack_from_test_data(env_tif):
    with rasterio.open(env_tif) as src:
        names = list(src.descriptions)
        data = src.read()
    # bio_2 normalized to bio2; numeric order
    assert names == [f"bio{i}" for i in range(1, 7)]
    assert data.shape == (6, 30, 40)
    # the NODATA corner (lon < 10.5, lat > 2.5) is NA in every band, the rest complete
    na = np.isnan(data).any(axis=0)
    assert na.sum() == 5 * 5
    assert na[:5, :5].all()
    assert (np.isnan(data).any(axis=0) == np.isnan(data).all(axis=0)).all()
    layers = (env_tif.parent / "layers.txt").read_text().split()
    assert layers == names


def test_env_stack_aggregate(tmp_path, env_tif):
    run_script("prep_env_stack.py", "--env-dir", os.path.join(TEST_DATA, "env"), "--pattern", "*.asc",
               "--clip", env_tif.parent / "area.geojson", "--aggregate", 2,
               "--out-tif", tmp_path / "agg.tif", "--out-layers", tmp_path / "layers.txt")
    with rasterio.open(tmp_path / "agg.tif") as src:
        assert (src.width, src.height) == (20, 15)
        assert src.res == pytest.approx((0.2, 0.2))


def test_env_stack_fails_without_layers(tmp_path, env_tif):
    p = run_script("prep_env_stack.py", "--env-dir", tmp_path, "--clip", env_tif.parent / "area.geojson",
                   "--out-tif", tmp_path / "x.tif", "--out-layers", tmp_path / "x.txt", check=False)
    assert p.returncode != 0 and "no files matching" in p.stderr


def test_extract_env_snaps_na_site(tmp_path, env_tif):
    sp = tmp_path / "sample_pop.tsv"
    sp.write_text("vcf_id\tpopid\nx1\t1\nx2\t1\nx3\t2\n")
    out = tmp_path / "env_sites.tsv"
    run_script("extract_env.py", "--latlongs", os.path.join(TEST_DATA, "spA", "spA_pop_latlongs.csv"),
               "--sample-pop", sp, "--env-tif", env_tif, "--max-snap-km", 10, "--out", out)
    e = pd.read_csv(out, sep="\t", dtype={"popid": str})
    assert e.popid.tolist() == ["1", "2"]  # only populations with samples
    assert e.n_ind.tolist() == [2, 1]
    assert 5 < e.snapped_km[0] < 10 and e.snapped_km[1] == 0
    assert not e.filter(like="bio").isna().any().any()


def test_extract_env_fails_when_snap_too_far(tmp_path, env_tif):
    sp = tmp_path / "sample_pop.tsv"
    sp.write_text("vcf_id\tpopid\nx1\t1\n")
    p = run_script("extract_env.py", "--latlongs", os.path.join(TEST_DATA, "spA", "spA_pop_latlongs.csv"),
                   "--sample-pop", sp, "--env-tif", env_tif, "--max-snap-km", 2, "--out", tmp_path / "o.tsv",
                   check=False)
    assert p.returncode != 0 and "no environmental data" in p.stderr


def test_get_boundaries_none_and_file(tmp_path, env_tif):
    area = env_tif.parent / "area.geojson"
    run_script("get_boundaries.py", "--source", "none", "--clip", area, "--cache-dir", tmp_path,
               "--out", tmp_path / "none.geojson")
    assert json.loads((tmp_path / "none.geojson").read_text())["features"] == []
    # a user-supplied boundary file is clipped to the study area
    src = tmp_path / "countries.geojson"
    src.write_text(json.dumps({"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"NAME": "Inside"},
         "geometry": {"type": "Polygon", "coordinates": [[[9, -1], [12, -1], [12, 4], [9, 4], [9, -1]]]}},
        {"type": "Feature", "properties": {"NAME": "Far away"},
         "geometry": {"type": "Polygon", "coordinates": [[[50, 50], [51, 50], [51, 51], [50, 50]]]}}]}))
    run_script("get_boundaries.py", "--source", src, "--clip", area, "--cache-dir", tmp_path,
               "--out", tmp_path / "b.geojson")
    names = [f["properties"]["name"] for f in json.loads((tmp_path / "b.geojson").read_text())["features"]]
    assert names == ["Inside"]
