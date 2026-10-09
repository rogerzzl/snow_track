"""Apparent lag and carry-over diagnostics of snowmelt-derived baseflow (Supplementary S9).

The routed snowmelt-derived surface runoff is used as a near-immediate melt
proxy and the routed snowmelt-derived subsurface runoff as the delayed
response.  A non-negative distributed-lag model with six lag bins is fitted by
non-negative least squares (needs ``scipy``); the bin coefficients are spread
uniformly into a daily kernel from which apparent mean / median / 90th
percentile lags and long-lag fractions follow.  These are model-derived
storage-release time scales, not isotope-equivalent water ages.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .metrics import annual_timing

AGE_BINS: List[Tuple[str, int, int]] = [
    ("0_30d", 0, 30),
    ("30_90d", 30, 90),
    ("90_180d", 90, 180),
    ("180_365d", 180, 365),
    ("1_2yr", 365, 730),
    ("2_5yr", 730, 1825),
]


def lag_bin_mean(x: np.ndarray, start: int, end: int) -> np.ndarray:
    """Mean of ``x`` over lags ``start <= tau < end`` for every time step (zero before the record)."""
    n = len(x)
    cs = np.concatenate([[0.0], np.cumsum(x)])
    t = np.arange(n)
    right = np.clip(t - start + 1, 0, n)
    left = np.clip(t - end + 1, 0, n)
    return (cs[right] - cs[left]) / float(end - start)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    order = np.argsort(values)
    v, w = values[order], weights[order]
    c = np.cumsum(w)
    if c[-1] <= 0:
        return np.nan
    return float(v[np.searchsorted(c, q * c[-1], side="left")])


def regression_metrics(obs: np.ndarray, pred: np.ndarray) -> Dict[str, float]:
    err = pred - obs
    ss_res = float(np.sum(err ** 2))
    ss_tot = float(np.sum((obs - obs.mean()) ** 2))
    return {
        "R2": 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan,
        "NSE": 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan,
        "RMSE": float(np.sqrt(np.mean(err ** 2))),
        "MAE": float(np.mean(np.abs(err))),
        "bias": float(np.mean(err)),
    }


def fit_distributed_lag(surf: np.ndarray, sub: np.ndarray, bins: Sequence[Tuple[str, int, int]] = AGE_BINS) -> Dict[str, object]:
    """Fit ``sub(t) ~ sum_j beta_j * mean_{tau in bin j} surf(t - tau)`` with beta >= 0."""
    try:
        from scipy.optimize import nnls
    except ImportError as exc:  # pragma: no cover
        raise ImportError("fit_distributed_lag needs scipy: pip install 'snowtrack[metrics]'") from exc
    surf = np.where(np.isfinite(surf), surf, 0.0).astype(float)
    y = np.asarray(sub, dtype=float)
    X = np.column_stack([lag_bin_mean(surf, s, e) for _, s, e in bins])
    valid = np.isfinite(y)
    beta, _ = nnls(X[valid], y[valid])
    fitted = X @ beta
    bsum = float(beta.sum())
    fractions = beta / bsum if bsum > 0 else np.full_like(beta, np.nan)

    kernel = pd.DataFrame([
        {"age_bin": lab, "lag_start_days": s, "lag_end_days": e, "beta_nnls": float(b), "apparent_fraction": float(f)}
        for (lab, s, e), b, f in zip(bins, beta, fractions)
    ])
    lags = np.concatenate([np.arange(s, e) for _, s, e in bins]).astype(float)
    h = np.concatenate([np.full(e - s, b / (e - s)) for (_, s, e), b in zip(bins, beta)])
    hsum = h.sum()
    if hsum > 0:
        p = h / hsum
        summary = {
            "apparent_mean_lag_days": float(np.sum(lags * p)),
            "apparent_median_lag_days": weighted_quantile(lags, p, 0.5),
            "apparent_tau90_days": weighted_quantile(lags, p, 0.9),
            "F_gt180d": float(p[lags >= 180].sum()),
            "F_gt365d": float(p[lags >= 365].sum()),
            "F_gt730d": float(p[lags >= 730].sum()),
        }
    else:
        summary = {k: np.nan for k in ("apparent_mean_lag_days", "apparent_median_lag_days", "apparent_tau90_days",
                                       "F_gt180d", "F_gt365d", "F_gt730d")}
    summary.update({f"model_{k}": v for k, v in regression_metrics(y[valid], fitted[valid]).items()})
    return {
        "beta": beta, "fractions": fractions, "kernel": kernel,
        "daily_kernel": pd.DataFrame({"lag_days": lags.astype(int), "h": h}),
        "fitted": fitted, "summary": summary,
    }


def allocate_origin_year(dates: pd.DatetimeIndex, surf: np.ndarray, sub: np.ndarray, daily_kernel: pd.DataFrame) -> pd.DataFrame:
    """Split each year's snowmelt-derived baseflow into same-year / previous-year / older apparent origins."""
    surf = np.where(np.isfinite(surf), surf, 0.0).astype(float)
    sub = np.where(np.isfinite(sub), sub, 0.0).astype(float)
    n = len(surf)
    year = dates.year.to_numpy()
    lags = daily_kernel["lag_days"].to_numpy(int)
    h = daily_kernel["h"].to_numpy(float)
    raw_total = np.zeros(n)
    same = np.zeros(n); prev = np.zeros(n); older = np.zeros(n); gt365 = np.zeros(n); gt730 = np.zeros(n)
    for lag, hh in zip(lags, h):
        if hh <= 0 or lag >= n:
            continue
        c = np.zeros(n)
        c[lag:] = hh * surf[: n - lag]
        raw_total += c
        oy = np.full(n, -1)
        oy[lag:] = year[: n - lag]
        dy = year - oy
        same += np.where(dy == 0, c, 0.0)
        prev += np.where(dy == 1, c, 0.0)
        older += np.where(dy >= 2, c, 0.0)
        if lag >= 365:
            gt365 += c
        if lag >= 730:
            gt730 += c
    scale = np.where(raw_total > 0, sub / np.where(raw_total > 0, raw_total, 1.0), 0.0)
    df = pd.DataFrame({
        "YEAR": year,
        "base_total_allocated": sub * (raw_total > 0),
        "base_same_year": same * scale, "base_previous_year": prev * scale, "base_older_years": older * scale,
        "base_gt365": gt365 * scale, "base_gt730": gt730 * scale,
    })
    g = df.groupby("YEAR").sum()
    den = g["base_total_allocated"].replace(0, np.nan)
    for k in ("same_year", "previous_year", "older_years", "gt365", "gt730"):
        g[f"F_{k}"] = g[f"base_{k}"] / den
    g["F_carryover"] = g["F_previous_year"] + g["F_older_years"]
    return g


def annual_subsurface_delay(std: pd.DataFrame, year_start_month: int = 4, start: Optional[int] = None,
                            end: Optional[int] = None, threshold: object = "p5") -> pd.DataFrame:
    """T50 and centroid delay of subsurface relative to surface snowmelt per melt-response year."""
    a_s = annual_timing(std["Q_SURF"], year_start_month, start, end, threshold)
    a_b = annual_timing(std["Q_SUB"], year_start_month, start, end, threshold)
    out = pd.DataFrame({
        "T50_surf": a_s["T50"], "T50_sub": a_b["T50"],
        "centroid_surf": a_s["CENTROID"], "centroid_sub": a_b["CENTROID"],
        "vol_surf": a_s["VOLUME"], "vol_sub": a_b["VOLUME"],
    })
    out["delay_T50_days"] = out["T50_sub"] - out["T50_surf"]
    out["delay_centroid_days"] = out["centroid_sub"] - out["centroid_surf"]
    return out


def next_winter_release(std: pd.DataFrame) -> pd.Series:
    """Ratio of subsurface snowmelt runoff in Jan-Mar of year y+1 to surface snowmelt runoff in Apr-Sep of year y."""
    y = std.index.year
    m = std.index.month
    surf = std["Q_SURF"][(m >= 4) & (m <= 9)].groupby(y[(m >= 4) & (m <= 9)]).sum()
    sub = std["Q_SUB"][(m >= 1) & (m <= 3)].groupby(y[(m >= 1) & (m <= 3)]).sum()
    sub.index = sub.index - 1
    ratio = sub / surf.replace(0, np.nan)
    ratio.name = "next_winter_release_ratio"
    return ratio.dropna()
