"""Canonical per-cell input schema and validation.

Every adapter converts model output into this schema; the tracker only ever
sees this schema.  One table = one model cell (grid cell, HRU, sub-basin or a
lumped catchment); one row = one model time step.

Time
----
``DATE`` (anything ``pandas.to_datetime`` understands) **or** ``YEAR, MONTH, DAY``
(optionally ``HOUR``).

Surface boundary (depth per time step, e.g. mm per day for daily output)
------------------------------------------------------------------------
``S_<SOURCE>``     water reaching the ground surface from each source
                   (``S_RAIN``, ``S_SNOW``, ``S_GLACIER``, ... any name).
``RUNOFF_SURF``    direct surface runoff leaving the cell.
``INFIL``          infiltration into the top soil layer (optional; default
                   ``max(sum(S_*) - RUNOFF_SURF, 0)``).

Soil column, for each layer n = 1..N (N detected from the ``W<n>_LIQ`` columns)
-------------------------------------------------------------------------------
``W<n>_LIQ``       liquid water storage at the END of the step (mm)      required
``W<n>_ICE``       soil-ice storage at the end of the step (mm)          optional
``Q<n>_PERC``      percolation out of the bottom of layer n (mm/step)    optional
``Q<n>_BASE``      drainage/baseflow from layer n to the channel         optional
``E<n>_TR``        evapotranspiration withdrawn from layer n             optional
``PHI<n>_FREEZE``  liquid -> ice phase change in layer n (mm/step)       optional
``PHI<n>_MELT``    ice -> liquid phase change in layer n (mm/step)       optional

If ``W<n>_ICE`` is present but the ``PHI`` columns are not, the net phase change
is inferred from the change in ice storage (``infer_phase_change=True``).
Percolation out of the bottom layer is treated as a loss from the column (deep
drainage that never reaches the channel).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

UNKNOWN = "UNKNOWN"
PREFERRED_SOURCE_ORDER: Tuple[str, ...] = ("RAIN", "SNOW", "GLACIER")
_LAYER_RE = re.compile(r"^W(\d+)_LIQ$")
_SOURCE_RE = re.compile(r"^S_([A-Z0-9]+)$")
_TIME_COLUMNS = ("YEAR", "MONTH", "DAY")


class SchemaError(ValueError):
    """Raised when an input table does not satisfy the canonical schema."""


def normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy whose column names are stripped and upper-cased."""
    out = df.copy()
    out.columns = [str(c).strip().upper() for c in out.columns]
    return out


def detect_sources(columns: Iterable[str]) -> List[str]:
    """Source names from ``S_<NAME>`` columns, rain/snow/glacier first."""
    found = []
    for c in columns:
        m = _SOURCE_RE.match(c)
        if m and m.group(1) not in ("TOT", "TOTAL", UNKNOWN):
            found.append(m.group(1))
    preferred = [s for s in PREFERRED_SOURCE_ORDER if s in found]
    others = sorted(s for s in found if s not in PREFERRED_SOURCE_ORDER)
    return preferred + others


def detect_layers(columns: Iterable[str]) -> int:
    """Number of soil layers N from ``W1_LIQ .. WN_LIQ``."""
    idx = sorted(int(_LAYER_RE.match(c).group(1)) for c in columns if _LAYER_RE.match(c))
    if not idx:
        raise SchemaError(
            "No soil-layer storage column found. At least 'W1_LIQ' is required "
            "(liquid soil-water storage at the end of each time step)."
        )
    if idx != list(range(1, len(idx) + 1)):
        raise SchemaError(f"Soil layers must be numbered 1..N without gaps; found W{idx}_LIQ.")
    return len(idx)


def parse_dates(df: pd.DataFrame) -> pd.DatetimeIndex:
    """Time stamps from ``DATE`` or ``YEAR/MONTH/DAY[/HOUR]``."""
    if "DATE" in df.columns:
        dates = pd.to_datetime(df["DATE"])
    elif all(c in df.columns for c in _TIME_COLUMNS):
        parts = {c.lower(): pd.to_numeric(df[c], errors="raise").astype(int) for c in _TIME_COLUMNS}
        if "HOUR" in df.columns:
            parts["hour"] = pd.to_numeric(df["HOUR"], errors="raise").astype(int)
        dates = pd.to_datetime(pd.DataFrame(parts))
    else:
        raise SchemaError("Need a DATE column or YEAR, MONTH, DAY columns.")
    dates = pd.DatetimeIndex(dates)
    if dates.has_duplicates:
        raise SchemaError("Duplicate time stamps found; one row per time step is required.")
    if not dates.is_monotonic_increasing:
        raise SchemaError("Rows must be sorted by time.")
    return dates


@dataclass
class CellData:
    """One model cell in the canonical schema, held as dense numpy arrays."""

    dates: pd.DatetimeIndex
    sources: Tuple[str, ...]
    n_layers: int
    supply: np.ndarray        # (T, K)  water reaching the ground, by source
    infil: np.ndarray         # (T,)
    runoff_surf: np.ndarray   # (T,)
    w_liq: np.ndarray         # (T, N)  end-of-step liquid storage
    w_ice: np.ndarray         # (T, N)  end-of-step ice storage
    q_perc: np.ndarray        # (T, N)
    q_base: np.ndarray        # (T, N)
    e_tr: np.ndarray          # (T, N)
    phi_freeze: np.ndarray    # (T, N)
    phi_melt: np.ndarray      # (T, N)
    has_ice: bool = False
    provided: Dict[str, object] = field(default_factory=dict)
    name: str = ""
    extra: Optional[pd.DataFrame] = None   # pass-through columns, not used by the tracker

    # ------------------------------------------------------------------ build
    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        name: str = "",
        infer_phase_change: bool = True,
        fill_missing: Optional[float] = None,
        keep_columns: Sequence[str] = (),
    ) -> "CellData":
        """Validate a canonical table and pack it into arrays."""
        df = normalise_columns(df).reset_index(drop=True)
        dates = parse_dates(df)
        sources = detect_sources(df.columns)
        if not sources:
            raise SchemaError(
                "No source supply column found. Provide at least one 'S_<SOURCE>' column, "
                "e.g. S_RAIN and S_SNOW (water reaching the ground from each source)."
            )
        n_layers = detect_layers(df.columns)
        if "RUNOFF_SURF" not in df.columns:
            raise SchemaError("Column 'RUNOFF_SURF' (direct surface runoff) is required.")

        T = len(df)
        provided: Dict[str, object] = {}

        def col(cname: str, default: Optional[float] = 0.0) -> np.ndarray:
            if cname in df.columns:
                provided[cname] = True
                values = pd.to_numeric(df[cname], errors="coerce").to_numpy(dtype=float)
                if np.isnan(values).any():
                    if fill_missing is None:
                        bad = int(np.isnan(values).sum())
                        raise SchemaError(
                            f"Column {cname} has {bad} missing/non-numeric values in {name or 'input'}. "
                            "Fix the input or pass fill_missing=0 (CLI: --fill-missing 0)."
                        )
                    values = np.where(np.isnan(values), fill_missing, values)
                return values
            provided[cname] = False
            if default is None:
                raise SchemaError(f"Column {cname} is required.")
            return np.full(T, default, dtype=float)

        supply = np.column_stack([col(f"S_{s}", None) for s in sources])
        runoff_surf = col("RUNOFF_SURF", None)
        if "INFIL" in df.columns:
            infil = col("INFIL")
        else:
            provided["INFIL"] = False
            infil = np.maximum(supply.sum(axis=1) - runoff_surf, 0.0)

        layers = range(1, n_layers + 1)
        w_liq = np.column_stack([col(f"W{n}_LIQ", None) for n in layers])
        w_ice = np.column_stack([col(f"W{n}_ICE") for n in layers])
        q_perc = np.column_stack([col(f"Q{n}_PERC") for n in layers])
        q_base = np.column_stack([col(f"Q{n}_BASE") for n in layers])
        e_tr = np.column_stack([col(f"E{n}_TR") for n in layers])
        phi_f = np.column_stack([col(f"PHI{n}_FREEZE") for n in layers])
        phi_m = np.column_stack([col(f"PHI{n}_MELT") for n in layers])
        has_ice = any(provided.get(f"W{n}_ICE", False) for n in layers)

        if has_ice and infer_phase_change:
            for j in range(n_layers):
                n = j + 1
                if provided.get(f"W{n}_ICE") and not (
                    provided.get(f"PHI{n}_FREEZE") or provided.get(f"PHI{n}_MELT")
                ):
                    net = np.zeros(T)
                    net[1:] = np.diff(w_ice[:, j])
                    phi_f[:, j] = np.maximum(net, 0.0)
                    phi_m[:, j] = np.maximum(-net, 0.0)
                    provided[f"PHI{n}_FREEZE"] = provided[f"PHI{n}_MELT"] = "inferred"

        if not np.any(q_base > 0):
            provided["_warning_no_baseflow"] = True

        extra = None
        keep = [c.upper() for c in keep_columns if c.upper() in df.columns]
        if keep:
            extra = df[keep].reset_index(drop=True)

        return cls(
            dates=dates, sources=tuple(sources), n_layers=n_layers,
            supply=supply, infil=infil, runoff_surf=runoff_surf,
            w_liq=w_liq, w_ice=w_ice, q_perc=q_perc, q_base=q_base, e_tr=e_tr,
            phi_freeze=phi_f, phi_melt=phi_m, has_ice=has_ice, provided=provided,
            name=name, extra=extra,
        )

    # ------------------------------------------------------------- utilities
    def __len__(self) -> int:
        return len(self.dates)

    @property
    def supply_total(self) -> np.ndarray:
        return self.supply.sum(axis=1)

    def slice(self, mask: np.ndarray, name: Optional[str] = None) -> "CellData":
        """Row subset (boolean mask), e.g. to drop input spin-up years."""
        return CellData(
            dates=self.dates[mask], sources=self.sources, n_layers=self.n_layers,
            supply=self.supply[mask], infil=self.infil[mask], runoff_surf=self.runoff_surf[mask],
            w_liq=self.w_liq[mask], w_ice=self.w_ice[mask], q_perc=self.q_perc[mask],
            q_base=self.q_base[mask], e_tr=self.e_tr[mask], phi_freeze=self.phi_freeze[mask],
            phi_melt=self.phi_melt[mask], has_ice=self.has_ice, provided=dict(self.provided),
            name=self.name if name is None else name,
            extra=None if self.extra is None else self.extra.loc[mask].reset_index(drop=True),
        )

    def budget_residuals(self) -> pd.DataFrame:
        """Per-step water-balance residual of the *bulk* inputs, layer by layer.

        ``residual_n(t) = [W_n(t) - W_n(t-1)] - [inflow_n - Q_n^perc - Q_n^base - E_n^tr]``
        with ``inflow_1 = INFIL`` and ``inflow_n = Q_{n-1}^perc``, where
        ``W_n = liquid + ice``.  A model whose exported fluxes close its own
        water balance gives ~0.  Non-zero values quantify fluxes the model did
        not export; the tracker absorbs them by rescaling the tracked storages to
        the bulk storage and, when needed, by assigning unexplained water to the
        UNKNOWN source (see ``docs/method.md``).
        """
        T, N = self.w_liq.shape
        total = self.w_liq + self.w_ice
        out = {"DATE": self.dates}
        out["SURFACE_RESIDUAL"] = self.supply_total - self.infil - self.runoff_surf
        for j in range(N):
            inflow = self.infil if j == 0 else self.q_perc[:, j - 1]
            dW = np.zeros(T)
            dW[1:] = np.diff(total[:, j])
            flux = inflow - self.q_perc[:, j] - self.q_base[:, j] - self.e_tr[:, j]
            res = dW - flux
            res[0] = np.nan
            out[f"LAYER{j + 1}_RESIDUAL"] = res
        return pd.DataFrame(out)

    def summary(self) -> Dict[str, object]:
        res = self.budget_residuals()
        return {
            "name": self.name,
            "steps": len(self),
            "start": str(self.dates[0].date()),
            "end": str(self.dates[-1].date()),
            "sources": list(self.sources),
            "n_layers": self.n_layers,
            "has_ice": self.has_ice,
            "provided": {k: v for k, v in self.provided.items() if not k.startswith("_")},
            "mean_abs_surface_residual": float(np.nanmean(np.abs(res["SURFACE_RESIDUAL"]))),
            "mean_abs_layer_residual": [
                float(np.nanmean(np.abs(res[f"LAYER{n}_RESIDUAL"]))) for n in range(1, self.n_layers + 1)
            ],
        }


def canonical_column_summary(
    n_layers: int = 3, sources: Sequence[str] = ("RAIN", "SNOW", "GLACIER")
) -> pd.DataFrame:
    """Table describing the canonical columns (used by ``snowtrack schema``)."""
    rows = [
        ("DATE or YEAR,MONTH,DAY[,HOUR]", "time stamp of the step", "required"),
        *[(f"S_{s}", f"{s.lower()} water reaching the ground surface (mm/step)", "at least one S_* column")
          for s in sources],
        ("RUNOFF_SURF", "direct surface runoff leaving the cell (mm/step)", "required"),
        ("INFIL", "infiltration into layer 1 (mm/step)", "optional: sum(S_*) - RUNOFF_SURF"),
    ]
    for n in range(1, n_layers + 1):
        rows += [
            (f"W{n}_LIQ", f"liquid storage of layer {n} at end of step (mm)", "required"),
            (f"W{n}_ICE", f"ice storage of layer {n} at end of step (mm)", "optional: 0"),
            (f"Q{n}_PERC", f"percolation out of layer {n} (mm/step)", "optional: 0"),
            (f"Q{n}_BASE", f"baseflow/drainage from layer {n} to channel (mm/step)", "optional: 0"),
            (f"E{n}_TR", f"evapotranspiration taken from layer {n} (mm/step)", "optional: 0"),
            (f"PHI{n}_FREEZE", f"liquid->ice in layer {n} (mm/step)", "optional: inferred from W_ICE"),
            (f"PHI{n}_MELT", f"ice->liquid in layer {n} (mm/step)", "optional: inferred from W_ICE"),
        ]
    return pd.DataFrame(rows, columns=["column", "meaning", "status"])
