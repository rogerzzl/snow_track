"""Aggregation of tracked per-cell outputs to annual / monthly / overall source fractions."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np
import pandas as pd

_SRC_RE = re.compile(r"^RUNOFF_([A-Z0-9]+)_TOTAL$")


def detect_output_sources(columns: Iterable[str]) -> List[str]:
    return [m.group(1) for c in columns if (m := _SRC_RE.match(c))]


def _safe_div(num, den):
    den = np.asarray(den, dtype=float)
    return np.where(np.abs(den) < 1e-12, np.nan, np.asarray(num, dtype=float) / np.where(np.abs(den) < 1e-12, 1.0, den))


def _sum_table(df: pd.DataFrame, keys: Sequence[str], sources: Sequence[str]) -> pd.DataFrame:
    cols = ["TOTAL_RUNOFF"] + [f"RUNOFF_{s}_{p}" for s in sources for p in ("TOTAL", "SURFACE", "SUBSURFACE")]
    g = df.groupby(list(keys), dropna=False)[cols].sum().reset_index()
    for s in sources:
        for p in ("TOTAL", "SURFACE", "SUBSURFACE"):
            g[f"F_{s}_{p}"] = _safe_div(g[f"RUNOFF_{s}_{p}"], g["TOTAL_RUNOFF"])
    return g


def aggregate(frames: Iterable[pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """Annual, monthly, overall and domain-climatology fraction tables.

    ``frames`` are tracker outputs, each with a ``GRID_ID`` column.
    """
    df = pd.concat(list(frames), ignore_index=True)
    if "GRID_ID" not in df.columns:
        df["GRID_ID"] = "cell"
    sources = detect_output_sources(df.columns)
    annual = _sum_table(df, ["GRID_ID", "YEAR"], sources)
    monthly = _sum_table(df, ["GRID_ID", "YEAR", "MONTH"], sources)

    overall_rows = []
    total = df["TOTAL_RUNOFF"].sum()
    for s in sources:
        overall_rows.append({
            "SOURCE": s,
            "RUNOFF_SOURCE_TOTAL": df[f"RUNOFF_{s}_TOTAL"].sum(),
            "RUNOFF_SOURCE_SURFACE": df[f"RUNOFF_{s}_SURFACE"].sum(),
            "RUNOFF_SOURCE_SUBSURFACE": df[f"RUNOFF_{s}_SUBSURFACE"].sum(),
            "TOTAL_RUNOFF": total,
            "F_SOURCE_TOTAL": float(_safe_div(df[f"RUNOFF_{s}_TOTAL"].sum(), total)),
            "F_SOURCE_SURFACE": float(_safe_div(df[f"RUNOFF_{s}_SURFACE"].sum(), total)),
            "F_SOURCE_SUBSURFACE": float(_safe_div(df[f"RUNOFF_{s}_SUBSURFACE"].sum(), total)),
        })
    overall = pd.DataFrame(overall_rows)

    # domain-wide monthly climatology (all cells and years pooled)
    clim = _sum_table(df, ["MONTH"], sources)
    return {"annual": annual, "monthly": monthly, "overall": overall, "monthly_climatology": clim}


def write_summaries(tables: Dict[str, pd.DataFrame], output_dir) -> Dict[str, Path]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    names = {
        "annual": "annual_source_fractions.csv",
        "monthly": "monthly_source_fractions.csv",
        "overall": "overall_source_fractions.csv",
        "monthly_climatology": "domain_monthly_climatology.csv",
    }
    written = {}
    for key, fname in names.items():
        if key in tables:
            p = out / fname
            tables[key].to_csv(p, index=False)
            written[key] = p
    return written


def load_tracked_dir(output_dir, pattern: str = "tracked_*") -> List[pd.DataFrame]:
    """Read all ``tracked_*`` files of a directory, adding ``GRID_ID``."""
    paths = sorted(Path(output_dir).glob(pattern))
    paths = [p for p in paths if p.suffix in (".csv", ".gz")]
    if not paths:
        raise FileNotFoundError(f"no {pattern} files in {output_dir}")
    frames = []
    for p in paths:
        name = p.name
        for suf in (".gz", ".csv"):
            if name.endswith(suf):
                name = name[: -len(suf)]
        name = name[len("tracked_"):] if name.startswith("tracked_") else name
        frames.append(pd.read_csv(p, parse_dates=["DATE"], compression="infer").assign(GRID_ID=name))
    return frames
