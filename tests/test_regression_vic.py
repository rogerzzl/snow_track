"""Reproduce the archived per-cell outputs of Zhu et al. from the raw VIC-glacier flux files."""

import numpy as np
import pandas as pd
import pytest

from snowtrack import CellData, SourceTracker, TrackerConfig, legacy
from snowtrack.adapters import get_adapter, load_config
from snowtrack.pipeline import RunOptions, run

CELLS = ["fluxes_31.0500_97.1500", "fluxes_32.1500_94.9500"]


def _expected(vic_example, name):
    return pd.read_csv(vic_example / "expected" / f"tracked_{name}.csv.gz", parse_dates=["DATE"])


@pytest.mark.parametrize("name", CELLS)
def test_reference_transcription_matches_archive(vic_example, name):
    raw = legacy.skip_input_years(legacy.read_vic_fluxes(vic_example / "input" / f"{name}.gz"), 2)
    ref = legacy.track_legacy(raw)
    exp = _expected(vic_example, name)
    assert len(ref) == len(exp)
    for c in ("RUNOFF_SNOW_SURFACE", "RUNOFF_SNOW_BASE", "F_W_SNOW", "F_W_UNKNOWN", "TOTAL_RUNOFF"):
        np.testing.assert_allclose(ref[c].to_numpy(), exp[c].to_numpy(), atol=1e-12, err_msg=c)


@pytest.mark.parametrize("name", CELLS)
def test_general_tracker_legacy_mode_matches_archive(vic_example, name):
    cfg = load_config(vic_example / "vic_classic_zhu2026.yaml")
    df = get_adapter(None, cfg).load(vic_example / "input" / f"{name}.gz")
    cell = CellData.from_dataframe(df, name=name)
    cell = cell.slice(np.asarray(cell.dates.year >= cell.dates[0].year + 2))
    out = SourceTracker(TrackerConfig(composition="legacy", store_storage_columns=False)).run(cell)
    exp = _expected(vic_example, name)
    assert (out["DATE"].to_numpy() == exp["DATE"].to_numpy()).all()
    np.testing.assert_allclose(out["RUNOFF_SNOW_SURFACE"], exp["RUNOFF_SNOW_SURFACE"], atol=1e-12)
    np.testing.assert_allclose(out["RUNOFF_SNOW_SUBSURFACE"], exp["RUNOFF_SNOW_BASE"], atol=1e-12)
    np.testing.assert_allclose(out["FW_SNOW"], exp["F_W_SNOW"], atol=1e-12)
    np.testing.assert_allclose(out["FW_UNKNOWN"], exp["F_W_UNKNOWN"], atol=1e-12)
    np.testing.assert_allclose(out["TOTAL_RUNOFF"], exp["TOTAL_RUNOFF"], atol=1e-12)


def test_other_compositions_stay_close_to_archive(vic_example):
    """'start' and 'end' differ from the archived convention only by the timing of mixing (< 1 % per year)."""
    name = CELLS[0]
    cfg = load_config(vic_example / "vic_classic_zhu2026.yaml")
    df = get_adapter(None, cfg).load(vic_example / "input" / f"{name}.gz")
    cell = CellData.from_dataframe(df, name=name)
    cell = cell.slice(np.asarray(cell.dates.year >= cell.dates[0].year + 2))
    exp = _expected(vic_example, name)
    exp_annual = exp.groupby("YEAR")["RUNOFF_SNOW_BASE"].sum()
    for comp in ("start", "end"):
        out = SourceTracker(TrackerConfig(composition=comp, store_storage_columns=False)).run(cell)
        np.testing.assert_allclose(out["RUNOFF_SNOW_SURFACE"], exp["RUNOFF_SNOW_SURFACE"], atol=1e-12)
        ann = out.groupby("YEAR")["RUNOFF_SNOW_SUBSURFACE"].sum()
        rel = np.abs(ann - exp_annual) / exp_annual.replace(0, np.nan)
        assert np.nanmax(rel) < 0.02


def test_cli_pipeline_reproduces_archive(vic_example, tmp_path):
    cfg = load_config(vic_example / "vic_classic_zhu2026.yaml")
    paths = run(vic_example / "input", tmp_path, "vic-classic", cfg,
                TrackerConfig(composition="legacy"), RunOptions(skip_input_years=2, limit=1, columns="compact"),
                log=lambda *a: None)
    assert len(paths) == 1
    out = pd.read_csv(paths[0], parse_dates=["DATE"])
    name = paths[0].name.replace("tracked_", "").replace(".csv", "")
    exp = _expected(vic_example, name)
    np.testing.assert_allclose(out["RUNOFF_SNOW_SUBSURFACE"], exp["RUNOFF_SNOW_BASE"], atol=1e-12)
    assert (tmp_path / "overall_source_fractions.csv").exists()
    assert (tmp_path / "tracking_run.json").exists()
