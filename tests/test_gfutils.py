"""gfutils.py: sample sheets, ID normalization, distances and raster sampling."""
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from gfutils import EnvStack, haversine_km, normalize_ids, read_latlongs, read_popmap


@pytest.mark.parametrize("text", [
    "sampleid popid\nA1 1\nA2 2\n",
    "sampleid\tpopid\nA1\t1\nA2\t2\n",
    "SampleID,PopID\nA1,1\nA2,2\n",
    pytest.param("sampleid  popid\nA1        1\nA2        2\n", id="aligned-columns"),
    pytest.param("sampleid\tpopid\nA1 \t1\nA2\t 2\n", id="mixed-tabs-spaces"),
    pytest.param("sampleid popid \nA1 1 \nA2 2\n", id="trailing-spaces"),
    pytest.param("sampleid popid\r\nA1 1\r\nA2 2\r\n", id="crlf"),
    pytest.param("SampleID, PopID\nA1, 1\nA2 ,2\n", id="comma-with-spaces"),
    pytest.param("sampleid popid\nA1 1\n\nA2 2\n", id="blank-line"),
])
def test_read_popmap_delimiters(tmp_path, text):
    f = tmp_path / "popmap.txt"
    f.write_text(text)
    pm = read_popmap(f)
    assert list(pm.columns) == ["sampleid", "popid"]
    assert pm.sampleid.tolist() == ["A1", "A2"]
    assert pm.popid.tolist() == ["1", "2"]


@pytest.mark.parametrize("text", [
    "sampleid popid\nA1 1\nA2\n",        # missing popid
    "sampleid,popid\nA1,1\n,2\n",        # empty sampleid
])
def test_read_popmap_empty_fields(tmp_path, text):
    f = tmp_path / "popmap.txt"
    f.write_text(text)
    with pytest.raises(SystemExit, match="line"):
        read_popmap(f)


def test_read_popmap_keeps_ids_as_text(tmp_path):
    f = tmp_path / "popmap.txt"
    f.write_text("sampleid popid\n007 01\n")
    pm = read_popmap(f)
    assert pm.sampleid.tolist() == ["007"] and pm.popid.tolist() == ["01"]


def test_read_popmap_missing_column(tmp_path):
    f = tmp_path / "bad.txt"
    f.write_text("id pop\nA1 1\n")
    with pytest.raises(SystemExit, match="sampleid"):
        read_popmap(f)


def test_read_latlongs_aliases_and_default_sitename(tmp_path):
    f = tmp_path / "ll.csv"
    f.write_text("popid,latitude,lon\n01,1.5,10.2\n2,2.5,11.0\n")
    ll = read_latlongs(f)
    assert list(ll.columns) == ["popid", "sitename", "lat", "long"]
    assert ll.popid.tolist() == ["01", "2"]  # kept as text, leading zero intact
    assert ll.sitename.tolist() == ["01", "2"]
    assert ll.lat.tolist() == [1.5, 2.5]


def test_read_latlongs_missing_column(tmp_path):
    f = tmp_path / "ll.csv"
    f.write_text("popid,lat\n1,1.5\n")
    with pytest.raises(SystemExit, match="long"):
        read_latlongs(f)


def test_normalize_ids():
    assert normalize_ids(["MCZA_1", "MCZA_2", "3"], "^MCZA_") == ["1", "2", "3"]
    assert normalize_ids(["MCZA_1"], "") == ["MCZA_1"]


def test_haversine_one_degree_latitude():
    assert haversine_km(10, 0, 10, 1) == pytest.approx(111.2, abs=0.2)
    assert haversine_km(10, 0, 10, 0) == 0


@pytest.fixture
def small_stack(tmp_path):
    """2 bands on a 4 x 4 grid of 0.1-degree cells (lon 10-10.4, lat 0-0.4);
    cell (row 0, col 0) is NA"""
    data = np.arange(32, dtype="float32").reshape(2, 4, 4)
    data[:, 0, 0] = np.nan
    path = tmp_path / "s.tif"
    with rasterio.open(path, "w", driver="GTiff", width=4, height=4, count=2, dtype="float32",
                       crs="EPSG:4326", transform=from_origin(10.0, 0.4, 0.1, 0.1), nodata=np.nan) as dst:
        dst.write(data)
        dst.set_band_description(1, "bio1")
        dst.set_band_description(2, "bio2")
    return EnvStack(path)


def test_envstack_basic(small_stack):
    s = small_stack
    assert s.names == ["bio1", "bio2"]
    assert s.inside(10.05, 0.35) and not s.inside(9.9, 0.2) and not s.inside(10.2, 0.5)
    vals, snap, n = s.sample(10.25, 0.15)  # row 2, col 2
    assert vals.tolist() == [10.0, 26.0] and snap == 0 and n == 1


def test_envstack_na_cell_and_snapping(small_stack):
    s = small_stack
    vals, snap, n = s.sample(10.05, 0.35)  # the NA cell, no snapping allowed
    assert np.isnan(vals).all() and n == 0
    vals, snap, n = s.sample(10.05, 0.35, max_snap_km=20)  # nearest valid: one cell away (~11 km)
    assert n == 1 and 10 < snap < 12
    assert vals[0] in (1.0, 4.0)  # east or south neighbour
    vals, snap, n = s.sample(10.05, 0.35, max_snap_km=5)  # too far
    assert np.isnan(vals).all() and snap > 5


def test_envstack_outside_and_buffer(small_stack):
    s = small_stack
    vals, snap, n = s.sample(9.0, 0.2)
    assert np.isnan(vals).all() and n == 0
    vals, snap, n = s.sample(10.25, 0.15, buffer_km=12)  # centre cell + 4 neighbours
    assert n == 5
    assert vals[0] == pytest.approx(np.mean([10, 6, 14, 9, 11]))
