"""Timing and pathway metrics used in Zhu et al. (Methods; Supplementary S4, S7, S8).

All functions work on an outlet (or cell) daily table with a ``DATE`` index and
the standard component columns produced by :func:`standardise_outlet_columns`:

``Q_SURF``   snowmelt-derived surface runoff        (``RUNOFF_SNOW_SURFACE[_OUTLET]``)
``Q_SUB``    snowmelt-derived subsurface runoff     (``RUNOFF_SNOW_SUBSURFACE[_OUTLET]`` / ``..._BASE``)
``Q_SNOW``   their sum
``Q_TOTAL``  total runoff                           (``TOTAL_RUNOFF[_OUTLET]``)
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

_CANDIDATES = {
    "Q_SURF": ["RUNOFF_{S}_SURFACE_OUTLET", "RUNOFF_{S}_SURFACE"],
    "Q_SUB": ["RUNOFF_{S}_SUBSURFACE_OUTLET", "RUNOFF_{S}_BASE_OUTLET", "RUNOFF_{S}_SUBSURFACE", "RUNOFF_{S}_BASE"],
    "Q_SNOW": ["RUNOFF_{S}_TOTAL_OUTLET", "RUNOFF_{S}_OUTLET", "RUNOFF_{S}_TOTAL", "RUNOFF_{S}"],
    "Q_TOTAL": ["TOTAL_RUNOFF_OUTLET", "TOTAL_RUNOFF"],
}


def read_outlet_series(path) -> pd.DataFrame:
    """Read an outlet CSV (``#`` comment lines allowed, ``.gz`` allowed)."""
    df = pd.read_csv(path, comment="#", parse_dates=["DATE"], compression="infer")
    return df.set_index("DATE").sort_index()


def standardise_outlet_columns(df: pd.DataFrame, source: str = "SNOW") -> pd.DataFrame:
    """Pick the component columns of ``source`` under any supported naming."""
    out = pd.DataFrame(index=df.index)
    cols = {c.upper(): c for c in df.columns}
    for key, cands in _CANDIDATES.items():
        for cand in cands:
            name = cand.replace("{S}", source.upper())
            if name in cols:
                out[key] = pd.to_numeric(df[cols[name]], errors="coerce").to_numpy(float)
                break
    if "Q_SNOW" not in out and {"Q_SURF", "Q_SUB"} <= set(out.columns):
        out["Q_SNOW"] = out["Q_SURF"] + out["Q_SUB"]
    missing = [k for k in ("Q_SURF", "Q_SUB", "Q_TOTAL") if k not in out]
    if missing:
        raise KeyError(f"could not find columns for {missing} (source {source}) in {list(df.columns)}")
    if "DATE" in out.columns:
        out = out.set_index("DATE")
    return out


# ------------------------------------------------------------------ timing
def t50_and_centroid(values: np.ndarray, day_index: Optional[np.ndarray] = None) -> Tuple[float, float, float]:
    """First day on which the cumulative sum reaches 50 % of the total, the centroid, and the total."""
    q = np.asarray(values, dtype=float)
    q = np.where(np.isfinite(q), q, 0.0)
    vol = float(q.sum())
    if vol <= 0:
        return np.nan, np.nan, vol
    x = np.arange(1, len(q) + 1, dtype=float) if day_index is None else np.asarray(day_index, dtype=float)
    csum = np.cumsum(q)
    idx = int(np.searchsorted(csum, 0.5 * vol, side="left"))
    return float(x[min(idx, len(x) - 1)]), float(np.sum(x * q) / vol), vol


def melt_year(dates: pd.DatetimeIndex, start_month: int = 1) -> np.ndarray:
    """Year label of a (melt) year that starts on day 1 of ``start_month``."""
    y = dates.year.to_numpy().copy()
    if start_month > 1:
        y[dates.month.to_numpy() < start_month] -= 1
    return y


def annual_timing(
    series: pd.Series,
    year_start_month: int = 1,
    start: Optional[int] = None,
    end: Optional[int] = None,
    threshold: object = None,
) -> pd.DataFrame:
    """Per-(melt-)year T50, centroid (day of the year, 1 = first day) and volume.

    ``threshold``: years with a volume at or below it get NaN timing.  ``None``
    (default, used for the calendar-year T50 trends of Table S4) excludes only
    zero volumes; ``"p5"`` = 5th percentile of the annual volumes (used for the
    melt-year lag diagnostics of S9); a float = absolute volume.
    """
    s = series.dropna()
    years = melt_year(s.index, year_start_month)
    rows = []
    for y in np.unique(years):
        grp = s[years == y]
        first = pd.Timestamp(int(y), year_start_month, 1)
        day = (grp.index - first).days.to_numpy() + 1
        t50, cen, vol = t50_and_centroid(grp.to_numpy(), day)
        rows.append({"YEAR": int(y), "T50": t50, "CENTROID": cen, "VOLUME": vol, "N_DAYS": len(grp)})
    out = pd.DataFrame(rows).set_index("YEAR")
    if start is not None:
        out = out[out.index >= start]
    if end is not None:
        out = out[out.index <= end]
    if threshold == "p5":
        vols = out["VOLUME"].replace(0, np.nan).dropna()
        thr = float(np.nanpercentile(vols, 5)) if len(vols) else 0.0
    elif threshold is None:
        thr = 0.0
    else:
        thr = float(threshold)
    bad = out["VOLUME"] <= thr
    out.loc[bad, ["T50", "CENTROID"]] = np.nan
    return out


def ols_trend(y: pd.Series) -> Dict[str, float]:
    """OLS slope per year with p-value (scipy) and n; NaNs dropped."""
    y = y.dropna()
    if len(y) < 3:
        return {"slope_per_year": np.nan, "p_value": np.nan, "n": int(len(y))}
    x = y.index.to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y.to_numpy(float), 1)
    p = np.nan
    try:
        from scipy import stats
        p = float(stats.linregress(x, y.to_numpy(float)).pvalue)
    except Exception:  # scipy optional
        pass
    return {"slope_per_year": float(slope), "p_value": p, "n": int(len(y))}


def timing_trends(
    std: pd.DataFrame,
    start: Optional[int] = None,
    end: Optional[int] = None,
    year_start_month: int = 1,
    metric: str = "T50",
    threshold: object = None,
) -> pd.DataFrame:
    """Trend of the timing metric (days per decade, **advance positive**) for each component."""
    rows = []
    for key in ("Q_SURF", "Q_SUB", "Q_SNOW", "Q_TOTAL"):
        if key not in std:
            continue
        at = annual_timing(std[key], year_start_month, start, end, threshold)
        tr = ols_trend(at[metric])
        rows.append({
            "component": key,
            "trend_days_per_decade_advance": -10.0 * tr["slope_per_year"],
            "p_value": tr["p_value"],
            "n_years": tr["n"],
            "mean_timing_doy": float(at[metric].mean()),
        })
    return pd.DataFrame(rows).set_index("component")


def buffering_index(std: pd.DataFrame, start: Optional[int] = None, end: Optional[int] = None,
                    year_start_month: int = 1, metric: str = "T50", threshold: object = None) -> Dict[str, float]:
    """Timing buffering index B = trend_surf - trend_total (days per decade, advance positive)."""
    tr = timing_trends(std, start, end, year_start_month, metric, threshold)["trend_days_per_decade_advance"]
    return {
        "trend_surf": float(tr.get("Q_SURF", np.nan)),
        "trend_sub": float(tr.get("Q_SUB", np.nan)),
        "trend_total": float(tr.get("Q_TOTAL", np.nan)),
        "B": float(tr.get("Q_SURF", np.nan) - tr.get("Q_TOTAL", np.nan)),
    }


# --------------------------------------------------------------- fractions
def monthly_fraction_climatology(std: pd.DataFrame, start: Optional[int] = None, end: Optional[int] = None) -> pd.DataFrame:
    """Monthly mean contribution (%) of surface and subsurface snowmelt to total runoff (Figure 1 style).

    ``*_MEAN`` averages the monthly fractions over the years, ``*_POOLED`` is the
    ratio of the pooled monthly sums.
    """
    df = std.copy()
    if start is not None:
        df = df[df.index.year >= start]
    if end is not None:
        df = df[df.index.year <= end]
    m = df.groupby([df.index.year, df.index.month])[["Q_SURF", "Q_SUB", "Q_TOTAL"]].sum()
    m.index.names = ["YEAR", "MONTH"]
    m = m[m["Q_TOTAL"] > 0]
    frac = pd.DataFrame({
        "F_SURF": m["Q_SURF"] / m["Q_TOTAL"],
        "F_SUB": m["Q_SUB"] / m["Q_TOTAL"],
    })
    out = frac.groupby(level="MONTH").mean() * 100.0
    out.columns = ["F_SURF_MEAN_PCT", "F_SUB_MEAN_PCT"]
    pooled = m.groupby(level="MONTH").sum()
    out["F_SURF_POOLED_PCT"] = 100.0 * pooled["Q_SURF"] / pooled["Q_TOTAL"]
    out["F_SUB_POOLED_PCT"] = 100.0 * pooled["Q_SUB"] / pooled["Q_TOTAL"]
    return out


def annual_fractions(std: pd.DataFrame) -> pd.DataFrame:
    a = std.groupby(std.index.year)[["Q_SURF", "Q_SUB", "Q_SNOW", "Q_TOTAL"]].sum()
    a.index.name = "YEAR"
    for k in ("Q_SURF", "Q_SUB", "Q_SNOW"):
        a[f"F_{k[2:]}"] = a[k] / a["Q_TOTAL"].replace(0, np.nan)
    return a


def late_season_lowflow(
    std: pd.DataFrame,
    months: Tuple[int, int] = (7, 10),
    window: int = 30,
    early: Tuple[int, int] = (1968, 1988),
    late: Tuple[int, int] = (1998, 2018),
    min_season_total: float = 1.0,
) -> Dict[str, object]:
    """Late-season subsurface snowmelt fraction and 30-day low-flow level (S8).

    Returns per-year values and the early/late period means and their change
    (``dF_sub`` in percentage points, ``dQ_low`` in the units of the series).
    """
    df = std.copy()
    rows = []
    for y, grp in df.groupby(df.index.year):
        season = grp[(grp.index.month >= months[0]) & (grp.index.month <= months[1])]
        if len(season) == 0:
            continue
        q_tot = float(season["Q_TOTAL"].sum())
        f_sub = float(season["Q_SUB"].sum() / q_tot) if q_tot >= min_season_total else np.nan
        ma = season["Q_TOTAL"].rolling(window, min_periods=window).mean()
        q_low = float(ma.min()) if ma.notna().any() else np.nan
        rows.append({"YEAR": int(y), "F_SUB_SEASON": f_sub, "Q_LOW_SEASON": q_low, "Q_SEASON_TOTAL": q_tot})
    per_year = pd.DataFrame(rows).set_index("YEAR")

    def mean_in(lo, hi, col):
        sel = per_year[(per_year.index >= lo) & (per_year.index <= hi)][col]
        return float(sel.mean()) if sel.notna().any() else np.nan

    f_early, f_late = mean_in(*early, "F_SUB_SEASON"), mean_in(*late, "F_SUB_SEASON")
    q_early, q_late = mean_in(*early, "Q_LOW_SEASON"), mean_in(*late, "Q_LOW_SEASON")
    return {
        "per_year": per_year,
        "F_sub_early_pct": 100.0 * f_early, "F_sub_late_pct": 100.0 * f_late,
        "dF_sub_pct_points": 100.0 * (f_late - f_early),
        "Q_low_early": q_early, "Q_low_late": q_late, "dQ_low": q_late - q_early,
    }
