"""Adapter for VIC 4.x classic-driver ASCII flux files (``fluxes_<lat>_<lon>``).

The classic driver writes one header-less whitespace-delimited file per grid
cell whose columns follow the ``OUTVAR`` list of the global parameter file
(after the time columns).  Variables that are written per soil layer
(``OUT_SOIL_MOIST``, ``OUT_SOIL_LIQ``, ``OUT_SOIL_ICE``, ...) expand into
``NLAYER`` columns.

Default configuration = the setup of Zhu et al. (see
``configs/vic_classic_zhu2026.yaml``)::

    outvars: [OUT_PREC, OUT_EVAP, OUT_RUNOFF, OUT_BASEFLOW, OUT_SWE,
              OUT_SNOW_MELT, OUT_RAINF, OUT_SNOWF, OUT_SOIL_MOIST]
    nlayer: 3
    time_columns: 3

Canonical mapping::

    S_RAIN      = RAINF              rain reaching the ground
    S_SNOW      = SNOW_MELT          seasonal snowmelt leaving the snowpack
    S_GLACIER   = <glacier column>   only if your VIC variant writes one
    RUNOFF_SURF = RUNOFF             VIC direct runoff
    INFIL       = max(S_tot - RUNOFF, 0)        (tracker default)
    W1_LIQ      = sum(SOIL_LIQ_n)  or sum(SOIL_MOIST_n) when no liquid/ice split is written
    W1_ICE      = sum(SOIL_ICE_n)  (only if OUT_SOIL_ICE was written)
    Q1_BASE     = BASEFLOW
    E1_TR       = EVAP - sublimation,  sublimation from the SWE balance
                  SWE(t) = SWE(t-1) + SNOWF(t) - SNOW_MELT(t) - SUBL(t)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd

from .base import Adapter

PER_LAYER_VARS = (
    "OUT_SOIL_MOIST", "OUT_SOIL_LIQ", "OUT_SOIL_ICE", "OUT_SOIL_TEMP",
    "OUT_SMLIQFRAC", "OUT_SMFROZFRAC",
)
DEFAULT_OUTVARS = [
    "OUT_PREC", "OUT_EVAP", "OUT_RUNOFF", "OUT_BASEFLOW", "OUT_SWE",
    "OUT_SNOW_MELT", "OUT_RAINF", "OUT_SNOWF", "OUT_SOIL_MOIST",
]


def expand_outvars(outvars: List[str], nlayer: int, per_layer=PER_LAYER_VARS) -> List[str]:
    """``OUT_SOIL_MOIST`` -> ``SOIL_MOIST1..N``; other names lose the ``OUT_`` prefix."""
    cols: List[str] = []
    for v in outvars:
        v = v.strip().upper()
        base = v[4:] if v.startswith("OUT_") else v
        if v in per_layer:
            cols.extend(f"{base}{n}" for n in range(1, nlayer + 1))
        else:
            cols.append(base)
    return cols


class VicClassicAdapter(Adapter):
    name = "vic-classic"
    patterns = ("fluxes_*", "fluxes_*.gz")

    def __init__(self, config=None) -> None:
        super().__init__(config)
        c = self.config
        self.outvars = list(c.get("outvars", DEFAULT_OUTVARS))
        self.nlayer = int(c.get("nlayer", 3))
        self.time_columns = int(c.get("time_columns", 3))
        self.columns = expand_outvars(self.outvars, self.nlayer)
        self.sources: Dict[str, str] = {k.upper(): v.upper() for k, v in
                                        c.get("sources", {"S_RAIN": "RAINF", "S_SNOW": "SNOW_MELT"}).items()}
        self.runoff_surf = c.get("runoff_surf", "RUNOFF").upper()
        self.baseflow = c.get("baseflow", "BASEFLOW").upper()
        self.soil_et = c.get("soil_et", "evap_minus_sublimation")
        self.sublimation = c.get("sublimation", "swe_balance")
        self.floor_negative_soil_et = bool(c.get("floor_negative_soil_et", False))
        self.keep = [k.upper() for k in c.get("keep", [])]

    # ------------------------------------------------------------------ io
    def read_raw(self, path: Path) -> pd.DataFrame:
        time_names = ["YEAR", "MONTH", "DAY", "HOUR"][: self.time_columns]
        names = time_names + self.columns
        df = pd.read_csv(path, sep=r"\s+", header=None, compression="infer", comment="#")
        if df.shape[1] != len(names):
            raise ValueError(
                f"{path}: expected {len(names)} columns ({', '.join(names)}) but found {df.shape[1]}. "
                "Set 'outvars', 'nlayer' and 'time_columns' in the adapter configuration to match "
                "the OUTVAR section of your VIC global parameter file."
            )
        df.columns = names
        return df

    def load(self, payload: Any) -> pd.DataFrame:
        raw = self.read_raw(Path(payload))
        return self.to_canonical(raw)

    # ----------------------------------------------------------- transform
    def _layer_sum(self, raw: pd.DataFrame, base: str):
        cols = [f"{base}{n}" for n in range(1, self.nlayer + 1) if f"{base}{n}" in raw.columns]
        if len(cols) == 0:
            return None
        return raw[cols].sum(axis=1).to_numpy(float)

    def to_canonical(self, raw: pd.DataFrame) -> pd.DataFrame:
        out = pd.DataFrame({c: raw[c] for c in ["YEAR", "MONTH", "DAY", "HOUR"] if c in raw.columns})
        for canon, column in self.sources.items():
            if column not in raw.columns:
                raise KeyError(f"source column {column} (for {canon}) not in VIC output columns {list(raw.columns)}")
            out[canon] = raw[column].to_numpy(float)
        out["RUNOFF_SURF"] = raw[self.runoff_surf].to_numpy(float)

        liq = self._layer_sum(raw, "SOIL_LIQ")
        ice = self._layer_sum(raw, "SOIL_ICE")
        moist = self._layer_sum(raw, "SOIL_MOIST")
        if liq is not None:
            out["W1_LIQ"] = liq
        elif moist is not None and ice is not None:
            out["W1_LIQ"] = moist - ice
        elif moist is not None:
            out["W1_LIQ"] = moist
        else:
            raise KeyError("VIC output must contain OUT_SOIL_MOIST or OUT_SOIL_LIQ")
        if ice is not None:
            out["W1_ICE"] = ice

        out["Q1_BASE"] = raw[self.baseflow].to_numpy(float)
        out["E1_TR"] = self._soil_et(raw)
        for k in self.keep:
            if k in raw.columns:
                out[k] = raw[k]
        return out

    def _sublimation(self, raw: pd.DataFrame) -> np.ndarray:
        if isinstance(self.sublimation, str) and self.sublimation.lower() == "swe_balance":
            swe = raw["SWE"].to_numpy(float)
            snowf = raw["SNOWF"].to_numpy(float)
            melt = raw["SNOW_MELT"].to_numpy(float)
            subl = np.zeros(len(raw))
            subl[1:] = np.maximum(swe[:-1] + snowf[1:] - melt[1:] - swe[1:], 0.0)
            return subl
        if self.sublimation in (None, "none", "zero"):
            return np.zeros(len(raw))
        col = str(self.sublimation).upper().replace("OUT_", "")
        return raw[col].to_numpy(float)

    def _soil_et(self, raw: pd.DataFrame) -> np.ndarray:
        mode = str(self.soil_et)
        if mode.lower() == "evap_minus_sublimation":
            et = raw["EVAP"].to_numpy(float) - self._sublimation(raw)
        elif mode.lower() == "evap":
            et = raw["EVAP"].to_numpy(float)
        else:
            et = raw[mode.upper().replace("OUT_", "")].to_numpy(float)
        if self.floor_negative_soil_et:
            et = np.maximum(et, 0.0)
        return et
