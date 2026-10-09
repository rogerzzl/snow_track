import numpy as np
import pandas as pd
import pytest

from snowtrack import CellData, SourceTracker, TrackerConfig
from snowtrack.adapters import get_adapter
from snowtrack.adapters.vic_classic import expand_outvars
from snowtrack.synthetic import make_synthetic_cell


def test_expand_outvars():
    cols = expand_outvars(["OUT_PREC", "OUT_SOIL_MOIST", "OUT_SOIL_ICE"], 3)
    assert cols == ["PREC", "SOIL_MOIST1", "SOIL_MOIST2", "SOIL_MOIST3", "SOIL_ICE1", "SOIL_ICE2", "SOIL_ICE3"]


def test_vic_adapter_builds_single_store(vic_example):
    ad = get_adapter("vic-classic", {})
    df = ad.load(vic_example / "input" / "fluxes_31.0500_97.1500.gz")
    assert list(df.columns) == ["YEAR", "MONTH", "DAY", "S_RAIN", "S_SNOW", "RUNOFF_SURF", "W1_LIQ", "Q1_BASE", "E1_TR"]
    raw = ad.read_raw(vic_example / "input" / "fluxes_31.0500_97.1500.gz")
    np.testing.assert_allclose(df["W1_LIQ"], raw[["SOIL_MOIST1", "SOIL_MOIST2", "SOIL_MOIST3"]].sum(axis=1))
    # sublimation from the SWE balance is subtracted from EVAP
    subl = np.zeros(len(raw))
    subl[1:] = np.maximum(raw.SWE.values[:-1] + raw.SNOWF.values[1:] - raw.SNOW_MELT.values[1:] - raw.SWE.values[1:], 0)
    np.testing.assert_allclose(df["E1_TR"], raw["EVAP"] - subl)


def test_vic_adapter_column_count_error(tmp_path):
    p = tmp_path / "fluxes_1_2"
    p.write_text("2000 1 1 0 0 0 0\n")
    with pytest.raises(ValueError, match="expected 14 columns"):
        get_adapter("vic-classic", {}).load(p)


def test_mapping_adapter_expressions_and_scale(tmp_path):
    raw = pd.DataFrame({
        "date": pd.date_range("2001-01-01", periods=4),
        "rain": [1, 0, 2, 0], "melt": [0, 3, 1, 0], "qs": [0.2, 0.5, 0.3, 0],
        "sm_a": [10, 11, 12, 12], "sm_b": [20, 21, 22, 22], "qb": [0.1, 0.1, 0.1, 0.1], "et": [0.5, 0.4, 0.6, 0.5],
    })
    raw.to_csv(tmp_path / "c1.csv", index=False)
    cfg = {"adapter": "mapping", "time": {"date": "date"},
           "columns": {"S_RAIN": "rain", "S_SNOW": "melt", "RUNOFF_SURF": "qs", "W1_LIQ": "sm_a + sm_b", "Q1_BASE": "qb", "E1_TR": "et"},
           "scale": {"S_RAIN": 2.0}}
    ad = get_adapter(None, cfg)
    items = list(ad.iter_items(tmp_path))
    assert [n for n, _ in items] == ["c1"]
    df = ad.load(items[0][1])
    np.testing.assert_allclose(df["W1_LIQ"], [30, 32, 34, 34])
    np.testing.assert_allclose(df["S_RAIN"], [2, 0, 4, 0])
    out = SourceTracker(TrackerConfig()).run(CellData.from_dataframe(df))
    assert len(out) == 4


def test_netcdf_adapter_roundtrip(tmp_path):
    xr = pytest.importorskip("xarray")
    pytest.importorskip("h5netcdf")
    base = make_synthetic_cell(years=1, seed=1, glacier=False, ice=True)
    time = pd.to_datetime(base[["YEAR", "MONTH", "DAY"]])
    T = len(base)
    lat, lon = [30.05, 30.15], [96.05, 96.15, 96.25]
    def grid(col, factor=1.0):
        a = np.empty((T, 2, 3)); a[:] = base[col].to_numpy()[:, None, None] * factor
        a[:, 1, 2] = np.nan        # one ocean/no-data cell
        return a
    liq = np.stack([grid("W1_LIQ"), grid("W2_LIQ"), grid("W3_LIQ")], axis=1)   # (time, lev, lat, lon)
    ice = np.stack([grid("W1_ICE"), grid("W2_ICE"), grid("W3_ICE")], axis=1)
    ds = xr.Dataset({
        "RAIN": (("time", "lat", "lon"), grid("S_RAIN") / 86400.0),
        "QSNOMELT": (("time", "lat", "lon"), grid("S_SNOW") / 86400.0),
        "QOVER": (("time", "lat", "lon"), grid("RUNOFF_SURF") / 86400.0),
        "QDRAI": (("time", "lat", "lon"), grid("Q3_BASE") / 86400.0),
        "SOILLIQ": (("time", "levsoi", "lat", "lon"), liq),
        "SOILICE": (("time", "levsoi", "lat", "lon"), ice),
    }, coords={"time": time, "lat": lat, "lon": lon, "levsoi": [0, 1, 2]})
    p = tmp_path / "clm.nc"
    ds.to_netcdf(p, engine="h5netcdf")
    cfg = {"adapter": "netcdf", "engine": "h5netcdf", "time_dim": "time", "cell_dims": ["lat", "lon"],
           "variables": {"S_RAIN": {"var": "RAIN", "scale": 86400}, "S_SNOW": {"var": "QSNOMELT", "scale": 86400},
                         "RUNOFF_SURF": {"var": "QOVER", "scale": 86400}, "Q1_BASE": {"var": "QDRAI", "scale": 86400},
                         "W1_LIQ": {"var": "SOILLIQ", "sum_over": "levsoi"}, "W1_ICE": {"var": "SOILICE", "sum_over": "levsoi"}}}
    ad = get_adapter(None, cfg)
    items = list(ad.iter_items(p))
    assert len(items) == 5 and items[0][0] == "cell_30.0500_96.0500"
    df = ad.load(items[0][1])
    np.testing.assert_allclose(df["W1_LIQ"], base[["W1_LIQ", "W2_LIQ", "W3_LIQ"]].sum(axis=1))
    np.testing.assert_allclose(df["S_SNOW"], base["S_SNOW"], rtol=1e-10)
    cell = CellData.from_dataframe(df)
    assert cell.has_ice and cell.provided["PHI1_FREEZE"] == "inferred"
    out = SourceTracker(TrackerConfig()).run(cell)
    assert np.nanmax(np.abs(out["RUNOFF_SOURCE_CLOSURE_ERROR"])) < 1e-9
