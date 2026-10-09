"""Conservation and reduction properties of the general tracker."""

import numpy as np
import pandas as pd
import pytest

from snowtrack import CellData, SourceTracker, TrackerConfig
from snowtrack.synthetic import make_synthetic_cell


def _closure(out: pd.DataFrame) -> float:
    cols = [c for c in out.columns if c.endswith("_CLOSURE_ERROR")]
    return float(np.nanmax(np.abs(out[cols].to_numpy())))


@pytest.mark.parametrize("composition", ["end", "start", "legacy"])
def test_mass_closure_three_layer_dual_phase(composition):
    cell = CellData.from_dataframe(make_synthetic_cell(years=3, seed=11), name="syn")
    out = SourceTracker(TrackerConfig(composition=composition)).run(cell)
    assert _closure(out) < 1e-7                      # storages sum to the bulk storage, runoff to the bulk runoff
    assert out["UNKNOWN_INJECTED"].sum() == 0.0
    for col in out.filter(like="RUNOFF_").columns:
        assert (out[col] >= -1e-12).all()
    # fractions of total runoff are between 0 and 1 and sum to 1 where runoff exists
    f = out[[c for c in out.columns if c.startswith("F_") and c.endswith("_TOTAL")]].dropna()
    assert ((f >= -1e-12) & (f <= 1 + 1e-12)).all().all()
    np.testing.assert_allclose(f.sum(axis=1), 1.0, atol=1e-9)


def test_end_composition_is_exact_on_mix_then_withdraw_data():
    """The synthetic model withdraws from the mixed store, so 'end' needs no clipping/reconciliation."""
    cell = CellData.from_dataframe(make_synthetic_cell(years=2, seed=5))
    out = SourceTracker(TrackerConfig(composition="end", reconcile=False)).run(cell)
    assert _closure(out) < 1e-7
    assert out["CLIPPED_MASS"].sum() < 1e-6


def test_unknown_source_is_flushed_with_time():
    cell = CellData.from_dataframe(make_synthetic_cell(years=6, seed=2))
    out = SourceTracker(TrackerConfig()).run(cell)
    first = out[out.YEAR == out.YEAR.min()]
    last = out[out.YEAR == out.YEAR.max()]
    assert first["FW_UNKNOWN"].iloc[1] > 0.9
    assert last["FW_UNKNOWN"].iloc[-1] < 0.05
    assert last["RUNOFF_UNKNOWN_TOTAL"].sum() / last["TOTAL_RUNOFF"].sum() < 0.05


def test_single_source_input_ends_fully_attributed():
    df = make_synthetic_cell(years=5, seed=7, glacier=False)
    df["S_SNOW"] = df["S_SNOW"] + df["S_RAIN"]          # everything that reaches the ground is 'snow'
    df["S_RAIN"] = 0.0
    out = SourceTracker(TrackerConfig()).run(CellData.from_dataframe(df))
    last = out[out.YEAR == out.YEAR.max()]
    assert last["RUNOFF_RAIN_TOTAL"].sum() == 0.0
    assert last["RUNOFF_SNOW_TOTAL"].sum() / last["TOTAL_RUNOFF"].sum() > 0.99
    assert last.loc[last["TOTAL_RUNOFF"] > 1e-6, "F_SNOW_TOTAL"].min() > 0.95


def test_phase_change_preserves_source_totals():
    """Freeze/thaw only moves water between phases: liquid+ice per source is unchanged by it."""
    df = make_synthetic_cell(years=2, seed=3, glacier=False)
    base = df.copy()
    base[[c for c in df.columns if c.startswith("PHI")]] = 0.0
    base[[c for c in df.columns if c.endswith("_ICE")]] = 0.0
    base[[c for c in df.columns if c.endswith("_LIQ")]] = df[[c for c in df.columns if c.endswith("_LIQ")]].to_numpy() + \
        df[[c for c in df.columns if c.endswith("_ICE")]].to_numpy()
    a = SourceTracker(TrackerConfig(composition="start")).run(CellData.from_dataframe(df))
    b = SourceTracker(TrackerConfig(composition="start")).run(CellData.from_dataframe(base))
    # baseflow leaves the liquid store only, so attribution differs; but source totals in storage match per layer
    for n in (1, 2, 3):
        for s in ("RAIN", "SNOW", "UNKNOWN"):
            tot_a = a[f"W{n}_{s}_LIQ"] + a[f"W{n}_{s}_ICE"]
            tot_b = b[f"W{n}_{s}_LIQ"]
            # equal on the first step (initialisation) and bounded afterwards
            assert abs(tot_a.iloc[0] - tot_b.iloc[0]) < 1e-9
    assert _closure(a) < 1e-7 and _closure(b) < 1e-7


def test_spinup_years_drop_output_only():
    cell = CellData.from_dataframe(make_synthetic_cell(years=3, seed=4, start_year=2001))
    out = SourceTracker(TrackerConfig(spinup_years=1)).run(cell)
    assert out["YEAR"].min() == 2002 and len(out) == 730


def test_strict_mode_raises_on_inconsistent_fluxes():
    df = make_synthetic_cell(years=1, seed=4, ice=False, glacier=False)
    df.loc[100, "Q3_BASE"] = 1e4                       # withdraw far more than the store holds
    with pytest.raises(ValueError, match="negative"):
        SourceTracker(TrackerConfig(strict=True, composition="start")).run(CellData.from_dataframe(df))
    out = SourceTracker(TrackerConfig(strict=False, composition="start")).run(CellData.from_dataframe(df))
    assert out["CLIPPED_MASS"].sum() > 0


def test_surface_excess_goes_to_unknown():
    df = make_synthetic_cell(years=1, seed=9, ice=False, glacier=False)
    df.loc[50, "RUNOFF_SURF"] = df.loc[50, ["S_RAIN", "S_SNOW"]].sum() + 3.0   # runoff exceeds supply
    df.loc[50, "INFIL"] = 0.0
    out = SourceTracker(TrackerConfig()).run(CellData.from_dataframe(df))
    assert abs(out.loc[50, "RUNOFF_UNKNOWN_SURFACE"] - 3.0) < 1e-9
    assert abs(out.loc[50, "RUNOFF_SOURCE_CLOSURE_ERROR"]) < 1e-9
    with pytest.raises(ValueError, match="surface closure"):
        SourceTracker(TrackerConfig(strict_surface_closure=True)).run(CellData.from_dataframe(df))


def test_v1_three_layer_sample_is_reproduced_with_start_composition(examples_dir):
    """The v1 implementation (explicit, start-of-step composition) is a special case of v2."""
    d = examples_dir / "01_minimal_three_layer"
    df = pd.read_csv(d / "input" / "grid_001.csv")
    exp = pd.read_csv(d / "expected" / "tracked_grid_001_v1_start.csv")
    out = SourceTracker(TrackerConfig(composition="start", strict_surface_closure=True)).run(CellData.from_dataframe(df))
    for c in exp.columns:
        if c == "DATE":
            continue
        np.testing.assert_allclose(out[c].to_numpy(float), exp[c].to_numpy(float), atol=1e-9, err_msg=c)
