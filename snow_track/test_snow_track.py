#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Focused checks for the S3 source-tracking implementation."""

from pathlib import Path

import numpy as np

from snow_track import DualPhaseSnowmeltTracker, TrackerConfig


ROOT = Path(__file__).resolve().parent
SAMPLE_INPUT = ROOT / "sampledata" / "input"
TEST_OUTPUT = ROOT / "sampledata" / "test_output"


def test_sample_conserves_storage_and_runoff_sources() -> None:
    tracker = DualPhaseSnowmeltTracker(
        SAMPLE_INPUT,
        TEST_OUTPUT,
        config=TrackerConfig(strict_surface_closure=True),
    )
    result = tracker.process_grid_cell(SAMPLE_INPUT / "grid_001.csv")

    closure_columns = [column for column in result.columns if column.endswith("_CLOSURE_ERROR")]
    max_storage_error = result[closure_columns].abs().to_numpy().max()
    max_runoff_error = result["RUNOFF_SOURCE_CLOSURE_ERROR"].abs().max()

    assert max_storage_error < 1.0e-9
    assert max_runoff_error < 1.0e-9
    assert np.isclose(result.loc[1, "RUNOFF_SNOW_SURFACE"], 1.0)
    assert result["RUNOFF_SNOW_SUBSURFACE"].iloc[-1] > 0.0
    finite_fraction = result["F_SNOW_TOTAL"].dropna()
    assert finite_fraction.between(0.0, 1.0).all()


if __name__ == "__main__":
    test_sample_conserves_storage_and_runoff_sources()
    print("All tests passed.")
