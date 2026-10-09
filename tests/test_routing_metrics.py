import numpy as np
import pandas as pd
import pytest

from snowtrack.lagdiag import AGE_BINS, annual_subsurface_delay, fit_distributed_lag, lag_bin_mean
from snowtrack.metrics import (annual_timing, buffering_index, late_season_lowflow, monthly_fraction_climatology,
                               read_outlet_series, standardise_outlet_columns, t50_and_centroid)
from snowtrack.routing import (RoutingConfig, haversine_km, inverse_gaussian_cdf, lag_steps, lohmann_kernel,
                               parse_latlon, route)


def _cell(name, lat, lon, n=10, surf=None, sub=None):
    dates = pd.date_range("2000-01-01", periods=n)
    surf = np.zeros(n) if surf is None else np.asarray(surf, float)
    sub = np.zeros(n) if sub is None else np.asarray(sub, float)
    df = pd.DataFrame({
        "DATE": dates, "RUNOFF_SNOW_SURFACE": surf, "RUNOFF_SNOW_SUBSURFACE": sub,
        "RUNOFF_SURF": surf, "BASEFLOW": sub, "TOTAL_RUNOFF": surf + sub,
    }).set_index("DATE")
    return (f"tracked_{name}_{lat:.4f}_{lon:.4f}", df)


# ------------------------------------------------------------------ routing
def test_parse_and_distance():
    assert parse_latlon("tracked_fluxes_31.0500_97.1500.csv") == (31.05, 97.15)
    assert parse_latlon("cell_7") is None
    assert abs(haversine_km(31.05, 97.15, 31.05, 97.15)) < 1e-9
    assert 9.4 < haversine_km(31.05, 97.05, 31.05, 97.15) < 9.6
    assert lag_steps(0.0, 1.0) == 0
    assert lag_steps(86.4, 1.0) == 1           # 86.4 km at 1 m/s = exactly 1 day


def test_inverse_gaussian_cdf_properties():
    mu, lam = 5.0 * 86400, 1e4 * 86400
    assert inverse_gaussian_cdf(0.0, mu, lam) == 0.0
    assert abs(inverse_gaussian_cdf(100 * mu, mu, lam) - 1.0) < 1e-9
    vals = [inverse_gaussian_cdf(t * 86400, mu, lam) for t in range(1, 40)]
    assert all(b >= a for a, b in zip(vals, vals[1:]))


@pytest.mark.parametrize("distance_km,velocity", [(50.0, 1.0), (300.0, 2.5), (120.0, 0.1)])
def test_lohmann_kernel_mass_and_mean_travel_time(distance_km, velocity):
    k = lohmann_kernel(distance_km * 1000, velocity, 800.0)
    assert abs(k.sum() - 1.0) < 1e-9 and (k >= 0).all()
    mean_days = float(np.sum(np.arange(len(k)) * k)) + 0.5         # bins are [k, k+1) days
    expected = distance_km * 1000 / velocity / 86400
    assert abs(mean_days - expected) < max(0.6, 0.02 * expected)
    assert len(lohmann_kernel(0.0, velocity, 800.0)) == 1


def test_lohmann_routing_conserves_mass_and_delays_baseflow():
    n = 400
    pulse = np.zeros(n); pulse[0] = 3.0
    cells = [_cell("a", 31.05, 97.15, n=n, surf=pulse, sub=pulse),
             _cell("b", 31.05, 98.95, n=n, surf=pulse, sub=pulse)]            # ~172 km from the outlet
    out = route(cells, outlet=(31.05, 97.15), config=RoutingConfig(surface_velocity=1.0, base_velocity=0.1, diffusivity=800.0))
    s = out["RUNOFF_SNOW_SURFACE_OUTLET"].to_numpy(); b = out["RUNOFF_SNOW_SUBSURFACE_OUTLET"].to_numpy()
    assert abs(s.sum() - 6.0) < 1e-6 and abs(b.sum() - 6.0) < 1e-6
    assert s[0] >= 3.0 and b[0] >= 3.0                                       # the outlet cell arrives immediately
    t_s = np.sum(np.arange(n) * (s - np.r_[3.0, np.zeros(n - 1)])) / 3.0    # mean arrival of the distant cell
    t_b = np.sum(np.arange(n) * (b - np.r_[3.0, np.zeros(n - 1)])) / 3.0
    assert 1.4 < t_s < 2.6 and 18.0 < t_b < 22.0                             # 172 km at 1 m/s ~ 2 d, at 0.1 m/s ~ 20 d
    np.testing.assert_allclose(out["TOTAL_RUNOFF_OUTLET"], s + b)
    assert out.attrs["n_cells"] == 2


def test_lag_scheme_shifts_components_by_their_own_lag():
    pulse = np.zeros(10); pulse[0] = 1.0
    cells = [_cell("a", 31.05, 97.15, surf=pulse, sub=pulse), _cell("b", 31.05, 98.95, surf=pulse, sub=pulse)]
    out = route(cells, outlet=(31.05, 97.15), config=RoutingConfig(surface_velocity=1.0, base_velocity=0.1, scheme="lag"))
    s = out["RUNOFF_SNOW_SURFACE_OUTLET"].to_numpy(); b = out["RUNOFF_SNOW_SUBSURFACE_OUTLET"].to_numpy()
    assert s[0] == 1.0 and s[2] == 1.0 and s.sum() == 2.0     # second cell: ~172 km at 1 m/s = 2 days
    assert b[0] == 1.0 and b[1:].sum() == 0.0                 # second cell's baseflow arrives after the record ends


def test_velocities_are_required():
    pulse = np.zeros(5); pulse[0] = 1.0
    with pytest.raises(ValueError, match="calibration parameter"):
        route([_cell("a", 31.05, 97.25, n=5, surf=pulse, sub=pulse)], outlet=(31.05, 97.15), config=RoutingConfig())


def test_route_with_cells_table_lags_and_weights():
    pulse = np.zeros(5); pulse[0] = 2.0
    cells = [_cell("x", 0, 0, n=5, surf=pulse, sub=pulse), _cell("y", 0, 0, n=5, surf=pulse, sub=pulse)]
    table = pd.DataFrame({"cell_id": ["x_0.0000_0.0000", "y_0.0000_0.0000"], "lag_surface": [0, 2], "lag_base": [1, 3], "weight": [1.0, 3.0]})
    out = route(cells, config=RoutingConfig(weight="area"), cells_table=table)
    np.testing.assert_allclose(out["RUNOFF_SNOW_SURFACE_OUTLET"].to_numpy(), [2 * 1 / 4, 0, 2 * 3 / 4, 0, 0])


def test_route_with_cells_table_distances():
    n = 200
    pulse = np.zeros(n); pulse[0] = 1.0
    cells = [_cell("p", 0, 0, n=n, surf=pulse, sub=pulse)]
    table = pd.DataFrame({"cell_id": ["p_0.0000_0.0000"], "distance_m": [86400.0 * 1.0 * 3]})   # 3 days at 1 m/s
    out = route(cells, config=RoutingConfig(surface_velocity=1.0, base_velocity=1.0, diffusivity=200.0), cells_table=table)
    s = out["RUNOFF_SNOW_SURFACE_OUTLET"].to_numpy()
    assert abs(s.sum() - 1.0) < 1e-6 and 2.0 < np.sum(np.arange(n) * s) < 3.6


# ------------------------------------------------------------------ metrics
def test_t50_and_annual_timing():
    q = np.zeros(10); q[3] = 1; q[4] = 1
    t50, cen, vol = t50_and_centroid(q)
    assert t50 == 4 and abs(cen - 4.5) < 1e-12 and vol == 2
    s = pd.Series(np.r_[np.zeros(100), np.ones(10), np.zeros(255)], index=pd.date_range("2001-01-01", periods=365))
    at = annual_timing(s, threshold=None)
    assert at.loc[2001, "T50"] == 105
    at4 = annual_timing(s, year_start_month=4, threshold=None)
    assert 2000 in at4.index and at4.loc[2001, "T50"] == 15


def test_lag_bin_mean_matches_naive():
    x = np.random.default_rng(0).random(50)
    out = lag_bin_mean(x, 2, 5)
    naive = np.array([np.mean([x[t - tau] if t - tau >= 0 else 0.0 for tau in range(2, 5)]) for t in range(50)])
    np.testing.assert_allclose(out, naive)


def test_distributed_lag_recovers_known_kernel():
    pytest.importorskip("scipy")
    rng = np.random.default_rng(1)
    surf = rng.gamma(1.0, 1.0, 4000)
    sub = 0.5 * lag_bin_mean(surf, 0, 30) + 0.3 * lag_bin_mean(surf, 180, 365)
    fit = fit_distributed_lag(surf, sub)
    beta = dict(zip([b[0] for b in AGE_BINS], fit["beta"]))
    assert abs(beta["0_30d"] - 0.5) < 1e-6 and abs(beta["180_365d"] - 0.3) < 1e-6
    assert fit["summary"]["model_R2"] > 0.999


def test_changdu_timing_trends_match_table_s4(outlets_dir):
    std = standardise_outlet_columns(read_outlet_series(outlets_dir / "changdu_snowmelt.csv.gz"))
    bi = buffering_index(std, 1968, 2018)
    assert abs(bi["trend_surf"] - 8.73) < 0.01       # Table S4, Changdu: surface snowmelt T50 advance (days per decade)
    assert abs(bi["trend_sub"] - 1.53) < 0.01        # Table S4, Changdu: subsurface snowmelt
    assert np.isfinite(bi["B"]) and bi["B"] > 0


def test_outlet_metric_chain_runs(outlets_dir):
    pytest.importorskip("scipy")
    std = standardise_outlet_columns(read_outlet_series(outlets_dir / "changdu_snowmelt.csv.gz"))
    std = std[(std.index.year >= 1968) & (std.index.year <= 2018)]
    clim = monthly_fraction_climatology(std)
    assert clim["F_SURF_MEAN_PCT"].idxmax() == 5            # surface snowmelt peaks in May at Changdu
    ls = late_season_lowflow(std)
    assert np.isfinite(ls["dF_sub_pct_points"]) and np.isfinite(ls["dQ_low"])
    d = annual_subsurface_delay(std, 4, 1968, 2017, threshold="p5")
    assert 50 < d["delay_T50_days"].median() < 110
    fit = fit_distributed_lag(std["Q_SURF"].to_numpy(), std["Q_SUB"].to_numpy())
    assert 0 < fit["summary"]["F_gt365d"] < 1 and fit["summary"]["model_R2"] > 0.5
