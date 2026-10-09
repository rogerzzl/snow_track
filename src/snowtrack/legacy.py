"""Reference transcription of the tracker that produced the archived results of Zhu et al.

This module is intentionally small and literal.  It re-implements, step by
step, the tracking code that was run on the VIC classic ``fluxes_*`` files of
the study (see ``docs/reproduce_paper.md``).
It exists so that the general tracker (:mod:`snowtrack.tracker`) can be shown
to reproduce the archived per-cell outputs to machine precision
(``tests/test_regression_vic.py``).  New work should use the general tracker.

Input: the 14-column VIC classic daily output used in the study
``YEAR MONTH DAY PREC EVAP RUNOFF BASEFLOW SWE SNOW_MELT RAINF SNOWF SOIL_MOIST1 SOIL_MOIST2 SOIL_MOIST3``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

VIC_COLUMNS = [
    "YEAR", "MONTH", "DAY", "PREC", "EVAP", "RUNOFF", "BASEFLOW",
    "SWE", "SNOW_MELT", "RAINF", "SNOWF", "SOIL_MOIST1", "SOIL_MOIST2", "SOIL_MOIST3",
]


def read_vic_fluxes(path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=None, names=VIC_COLUMNS)
    df["DATE"] = pd.to_datetime(df[["YEAR", "MONTH", "DAY"]])
    return df


def skip_input_years(df: pd.DataFrame, years: int) -> pd.DataFrame:
    """Drop the first ``years`` calendar years *before* tracking (archived convention)."""
    if years <= 0 or df["YEAR"].nunique() <= years:
        return df.reset_index(drop=True)
    first = int(df["YEAR"].iloc[0]) + years
    return df.loc[df["DATE"] >= pd.Timestamp(first, 1, 1)].reset_index(drop=True)


def track_legacy(df: pd.DataFrame) -> pd.DataFrame:
    """Rain/snow/unknown tracker of the archived run (literal transcription)."""
    evap = df["EVAP"].to_numpy(float)
    runoff = df["RUNOFF"].to_numpy(float)
    basef = df["BASEFLOW"].to_numpy(float)
    swe = df["SWE"].to_numpy(float)
    melt = df["SNOW_MELT"].to_numpy(float)
    rainf = df["RAINF"].to_numpy(float)
    snowf = df["SNOWF"].to_numpy(float)
    W = (df["SOIL_MOIST1"] + df["SOIL_MOIST2"] + df["SOIL_MOIST3"]).to_numpy(float)
    n = len(df)

    subl = np.zeros(n)
    subl[1:] = np.maximum(swe[:-1] + snowf[1:] - melt[1:] - swe[1:], 0.0)
    infil = np.maximum(rainf + melt - runoff, 0.0)

    f_snow = np.zeros(n)
    f_rain = np.zeros(n)
    f_unk = np.ones(n)
    f_i = np.zeros(n)
    r_snow = np.zeros(n)
    b_snow = np.zeros(n)

    for t in range(1, n):
        if melt[t] > 0:
            f_i[t] = melt[t] / (melt[t] + rainf[t]) if rainf[t] > 0 else 1.0
            r_snow[t] = max(0.0, melt[t] - f_i[t] * infil[t])
        et = evap[t] - subl[t]
        if W[t] > 0:
            out = (et + basef[t])
            fs = (f_snow[t - 1] * W[t - 1] + f_i[t] * infil[t] - f_snow[t - 1] * out) / W[t]
            fr = (f_rain[t - 1] * W[t - 1] + (1.0 - f_i[t]) * infil[t] - f_rain[t - 1] * out) / W[t]
            fu = (f_unk[t - 1] * W[t - 1] - f_unk[t - 1] * out) / W[t]
            fs, fr, fu = (min(max(v, 0.0), 1.0) for v in (fs, fr, fu))
            tot = fs + fr + fu
            if tot > 0:
                fs, fr, fu = fs / tot, fr / tot, fu / tot
            f_snow[t], f_rain[t], f_unk[t] = fs, fr, fu
        b_snow[t] = f_snow[t] * basef[t]

    q_snow = r_snow + b_snow
    total = runoff + basef
    f_q = np.where(total > 0, q_snow / np.where(total > 0, total, 1.0), 0.0)
    return pd.DataFrame({
        "DATE": df["DATE"].to_numpy(),
        "YEAR": df["YEAR"].to_numpy(),
        "MONTH": df["MONTH"].to_numpy(),
        "DAY": df["DAY"].to_numpy(),
        "TOTAL_RUNOFF": total,
        "RUNOFF_SNOW": q_snow,
        "RUNOFF_SNOW_SURFACE": r_snow,
        "RUNOFF_SNOW_BASE": b_snow,
        "F_Q_SNOW": f_q,
        "F_W_SNOW": f_snow,
        "F_W_RAIN": f_rain,
        "F_W_UNKNOWN": f_unk,
    })
