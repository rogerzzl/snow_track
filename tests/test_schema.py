import numpy as np
import pandas as pd
import pytest

from snowtrack import CellData, SchemaError
from snowtrack.synthetic import make_synthetic_cell


def _minimal(n=5):
    return pd.DataFrame({
        "date": pd.date_range("2000-01-01", periods=n, freq="D"),
        "s_rain": [0, 2, 0, 1, 0], "S_SNOW": [0, 1, 3, 0, 0],
        "runoff_surf": [0, 1, 1, 0, 0], "W1_LIQ": [50, 51, 52, 52, 51], "q1_base": [0, 0.5, 0.5, 0.5, 0.5],
    })


def test_detects_sources_layers_and_defaults():
    cell = CellData.from_dataframe(_minimal(), name="m")
    assert cell.sources == ("RAIN", "SNOW")
    assert cell.n_layers == 1 and not cell.has_ice
    # INFIL defaults to sum(S) - RUNOFF_SURF, floored at zero
    np.testing.assert_allclose(cell.infil, [0, 2, 2, 1, 0])
    assert cell.provided["INFIL"] is False and cell.provided["E1_TR"] is False
    assert np.all(cell.e_tr == 0)


def test_three_layer_dual_phase_detection():
    df = make_synthetic_cell(years=1)
    cell = CellData.from_dataframe(df)
    assert cell.n_layers == 3 and cell.has_ice and cell.sources == ("RAIN", "SNOW", "GLACIER")
    res = cell.budget_residuals()
    assert np.nanmax(np.abs(res.iloc[:, 1:].to_numpy())) < 1e-9


def test_phase_change_inferred_from_ice_storage():
    df = _minimal()
    df["W1_ICE"] = [10, 12, 9, 9, 9]
    cell = CellData.from_dataframe(df)
    assert cell.provided["PHI1_FREEZE"] == "inferred"
    np.testing.assert_allclose(cell.phi_freeze[:, 0], [0, 2, 0, 0, 0])
    np.testing.assert_allclose(cell.phi_melt[:, 0], [0, 0, 3, 0, 0])


def test_errors_are_informative():
    with pytest.raises(SchemaError, match="S_<SOURCE>"):
        CellData.from_dataframe(_minimal().drop(columns=["s_rain", "S_SNOW"]))
    with pytest.raises(SchemaError, match="W1_LIQ"):
        CellData.from_dataframe(_minimal().drop(columns=["W1_LIQ"]))
    with pytest.raises(SchemaError, match="RUNOFF_SURF"):
        CellData.from_dataframe(_minimal().drop(columns=["runoff_surf"]))
    df = _minimal(); df.loc[2, "W1_LIQ"] = np.nan
    with pytest.raises(SchemaError, match="missing"):
        CellData.from_dataframe(df)
    assert CellData.from_dataframe(df, fill_missing=0.0).w_liq[2, 0] == 0.0


def test_layer_numbering_must_be_contiguous():
    df = _minimal().rename(columns={"W1_LIQ": "W2_LIQ"})
    with pytest.raises(SchemaError, match="numbered"):
        CellData.from_dataframe(df)
