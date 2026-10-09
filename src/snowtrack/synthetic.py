"""Synthetic, exactly closed test cells (used by the demo and the tests).

A toy cold-region bucket model: seasonal rain and snowfall, degree-day
snowmelt, a small summer glacier-melt term, a three-layer soil column with
percolation, baseflow from the bottom layer, layer-wise transpiration and a
freeze/thaw cycle of soil water.  The model is integrated in the same order the
tracker uses (transport, then phase change), so the exported storages close the
exported fluxes exactly and the tracker's closure errors are zero to machine
precision.  The output is a canonical-schema DataFrame.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def make_synthetic_cell(
    years: int = 4,
    start_year: int = 2001,
    seed: int = 0,
    n_layers: int = 3,
    ice: bool = True,
    glacier: bool = True,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    dates = pd.date_range(f"{start_year}-01-01", f"{start_year + years - 1}-12-31", freq="D")
    T = len(dates)
    doy = dates.dayofyear.to_numpy()
    # temperature and precipitation climatology (degrees C, mm/day)
    temp = -8.0 + 20.0 * np.sin(2 * np.pi * (doy - 110) / 365.25) + rng.normal(0, 2.5, T)
    precip = np.maximum(rng.gamma(0.6, 4.0, T) * (0.4 + 0.8 * np.clip(np.sin(2 * np.pi * (doy - 90) / 365.25), 0, 1)), 0)
    rain = np.where(temp > 1.0, precip, 0.0)
    snowfall = np.where(temp <= 1.0, precip, 0.0)

    cap = np.array([40.0, 120.0, 250.0])[:n_layers]
    k_perc = np.array([0.08, 0.03, 0.0])[:n_layers]
    k_base = np.zeros(n_layers); k_base[-1] = 0.012
    pet_frac = np.array([0.5, 0.35, 0.15])[:n_layers]

    swe = 0.0
    glac_ice = 1.0e6
    w_liq = cap * 0.5
    w_ice = np.zeros(n_layers)
    cols = {k: np.zeros(T) for k in ["S_RAIN", "S_SNOW", "S_GLACIER", "INFIL", "RUNOFF_SURF", "SWE"]}
    lay = {f"{p}{n}": np.zeros(T) for n in range(1, n_layers + 1)
           for p in ("W", "WI", "QP", "QB", "ET", "PF", "PM")}

    for t in range(T):
        # snowpack
        swe += snowfall[t]
        melt = min(swe, max(temp[t], 0.0) * 3.0) if t > 0 else 0.0
        swe -= melt
        gmelt = (max(temp[t] - 4.0, 0.0) * 1.0 * (doy[t] > 150) * (doy[t] < 280)) if glacier else 0.0
        gmelt = min(gmelt, glac_ice); glac_ice -= gmelt
        s_tot = rain[t] + melt + gmelt
        # surface partition (VIC-like saturation excess)
        sat = w_liq[0] / cap[0]
        r_surf = s_tot * (0.15 + 0.6 * sat ** 3) if t > 0 else 0.0
        infil = s_tot - r_surf
        # transport, top-down
        inflow = infil
        pet = max(0.3 + 3.5 * max(temp[t], 0.0) / 20.0, 0.0)
        for n in range(n_layers):
            if t == 0:
                break
            avail = w_liq[n] + inflow
            et = min(pet * pet_frac[n] * (avail / cap[n]), avail * 0.5)
            perc = min(k_perc[n] * avail * (avail / cap[n]), avail - et)
            base = min(k_base[n] * avail, avail - et - perc)
            w_liq[n] = avail - et - perc - base
            lay[f"QP{n + 1}"][t], lay[f"QB{n + 1}"][t], lay[f"ET{n + 1}"][t] = perc, base, et
            inflow = perc
        # phase change
        if ice and t > 0:
            for n in range(n_layers):
                depth_damp = 0.6 ** n
                if temp[t] < -1.0:
                    fr = min(0.04 * depth_damp * (-temp[t] - 1.0), 1.0) * w_liq[n]
                    me = 0.0
                else:
                    fr = 0.0
                    me = min(0.08 * depth_damp * max(temp[t], 0.0), 1.0) * w_ice[n]
                w_liq[n] += me - fr
                w_ice[n] += fr - me
                lay[f"PF{n + 1}"][t], lay[f"PM{n + 1}"][t] = fr, me
        cols["S_RAIN"][t], cols["S_SNOW"][t], cols["S_GLACIER"][t] = (rain[t], melt, gmelt) if t > 0 else (0, 0, 0)
        cols["INFIL"][t], cols["RUNOFF_SURF"][t], cols["SWE"][t] = infil, r_surf, swe
        for n in range(n_layers):
            lay[f"W{n + 1}"][t] = w_liq[n]
            lay[f"WI{n + 1}"][t] = w_ice[n]

    out = pd.DataFrame({"YEAR": dates.year, "MONTH": dates.month, "DAY": dates.day})
    out["S_RAIN"] = cols["S_RAIN"]
    out["S_SNOW"] = cols["S_SNOW"]
    if glacier:
        out["S_GLACIER"] = cols["S_GLACIER"]
    out["INFIL"] = cols["INFIL"]
    out["RUNOFF_SURF"] = cols["RUNOFF_SURF"]
    for n in range(1, n_layers + 1):
        out[f"W{n}_LIQ"] = lay[f"W{n}"]
        if ice:
            out[f"W{n}_ICE"] = lay[f"WI{n}"]
        out[f"Q{n}_PERC"] = lay[f"QP{n}"]
        out[f"Q{n}_BASE"] = lay[f"QB{n}"]
        out[f"E{n}_TR"] = lay[f"ET{n}"]
        if ice:
            out[f"PHI{n}_FREEZE"] = lay[f"PF{n}"]
            out[f"PHI{n}_MELT"] = lay[f"PM{n}"]
    out["SWE"] = cols["SWE"]
    return out
