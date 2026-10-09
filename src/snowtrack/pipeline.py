"""Run the tracker over many cells (directory of files or a netCDF dataset)."""

from __future__ import annotations

import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from . import __version__
from .adapters import get_adapter
from .aggregate import aggregate, write_summaries
from .schema import CellData
from .tracker import SourceTracker, TrackerConfig

COLUMN_PROFILES = ("full", "standard", "compact")


@dataclass
class RunOptions:
    skip_input_years: int = 0          # drop the first N years of *input* before tracking
    fill_missing: Optional[float] = None
    keep_columns: Sequence[str] = ()
    columns: str = "standard"          # full | standard | compact
    gzip: bool = False
    float_format: Optional[str] = None
    workers: int = 1
    summarize: bool = True
    limit: Optional[int] = None
    infer_phase_change: bool = True
    overwrite: bool = True

    def __post_init__(self) -> None:
        if self.columns not in COLUMN_PROFILES:
            raise ValueError(f"columns must be one of {COLUMN_PROFILES}")


def build_cell(df: pd.DataFrame, name: str, opts: RunOptions) -> CellData:
    cell = CellData.from_dataframe(
        df, name=name, infer_phase_change=opts.infer_phase_change,
        fill_missing=opts.fill_missing, keep_columns=opts.keep_columns,
    )
    if opts.skip_input_years > 0:
        first = int(cell.dates[0].year) + opts.skip_input_years
        mask = cell.dates.year >= first
        if mask.sum() == 0:
            raise ValueError(f"{name}: skip_input_years={opts.skip_input_years} leaves no data")
        cell = cell.slice(np.asarray(mask))
    return cell


def select_columns(out: pd.DataFrame, profile: str, sources: Sequence[str]) -> pd.DataFrame:
    if profile == "full":
        return out
    base = ["DATE", "YEAR", "MONTH", "DAY", "TOTAL_RUNOFF", "RUNOFF_SURF", "BASEFLOW"]
    runoff = [f"RUNOFF_{s}_{p}" for s in sources for p in ("TOTAL", "SURFACE", "SUBSURFACE")]
    if profile == "compact":
        keep = base + runoff
    else:
        keep = base + runoff
        keep += [f"F_{s}_{p}" for s in sources for p in ("TOTAL", "SURFACE", "SUBSURFACE")]
        keep += [f"FW_{s}" for s in sources]
        keep += [c for c in out.columns if c.startswith("BUDGET_RESIDUAL_")]
        keep += ["CLIPPED_MASS", "UNKNOWN_INJECTED", "RUNOFF_SOURCE_CLOSURE_ERROR"]
    extra = [c for c in out.columns if c not in keep and not _is_internal(c)]
    keep = [c for c in keep if c in out.columns] + extra
    return out[keep]


def _is_internal(c: str) -> bool:
    return (
        c.startswith(("W", "F_", "FW_", "BUDGET_RESIDUAL_", "RUNOFF_")) and ("_" in c)
        or c in ("CLIPPED_MASS", "UNKNOWN_INJECTED", "RUNOFF_SOURCE_CLOSURE_ERROR")
    )


def track_dataframe(df: pd.DataFrame, name: str, tracker_config: TrackerConfig, opts: RunOptions) -> pd.DataFrame:
    cell = build_cell(df, name, opts)
    cfg = tracker_config
    if opts.columns != "full" and cfg.store_storage_columns:
        cfg = TrackerConfig(**{**cfg.__dict__, "store_storage_columns": False})
    out = SourceTracker(cfg).run(cell)
    sources = list(cell.sources) + ["UNKNOWN"]
    return select_columns(out, opts.columns, sources)


def _worker(args: Tuple) -> Dict[str, Any]:
    name, payload, adapter_name, adapter_config, tracker_config, opts, output_dir = args
    t0 = time.time()
    adapter = get_adapter(adapter_name, adapter_config)
    df = adapter.load(payload)
    out = track_dataframe(df, name, tracker_config, opts)
    out_path = Path(output_dir) / f"tracked_{name}.csv{'.gz' if opts.gzip else ''}"
    out.to_csv(out_path, index=False, float_format=opts.float_format, compression="gzip" if opts.gzip else None)
    diag = {
        "cell": name,
        "rows": int(len(out)),
        "max_abs_closure_error": float(np.nanmax(np.abs(out["RUNOFF_SOURCE_CLOSURE_ERROR"]))) if "RUNOFF_SOURCE_CLOSURE_ERROR" in out else None,
        "unknown_fraction_last_year": _unknown_last_year(out),
        "seconds": round(time.time() - t0, 2),
    }
    return {"path": str(out_path), "diag": diag}


def _unknown_last_year(out: pd.DataFrame) -> Optional[float]:
    if "RUNOFF_UNKNOWN_TOTAL" not in out or len(out) == 0:
        return None
    last = out["YEAR"].iloc[-1]
    sel = out[out["YEAR"] == last]
    tot = sel["TOTAL_RUNOFF"].sum()
    return float(sel["RUNOFF_UNKNOWN_TOTAL"].sum() / tot) if tot > 0 else None


def run(
    source,
    output_dir,
    adapter_name: Optional[str] = None,
    adapter_config: Optional[dict] = None,
    tracker_config: Optional[TrackerConfig] = None,
    options: Optional[RunOptions] = None,
    log=print,
) -> List[Path]:
    """Track every cell of ``source`` and write ``tracked_<cell>.csv`` files plus summaries."""
    opts = options or RunOptions()
    tcfg = tracker_config or TrackerConfig()
    adapter_config = dict(adapter_config or {})
    adapter = get_adapter(adapter_name, adapter_config)
    adapter_name = adapter.name
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    items = list(adapter.iter_items(source))
    if opts.limit:
        items = items[: opts.limit]
    if not items:
        raise FileNotFoundError(f"{adapter.describe()}: no cells found in {source}")
    log(f"snowtrack {__version__}: {len(items)} cell(s) via {adapter.describe()}, "
        f"composition={tcfg.composition}, workers={opts.workers}")

    jobs = [(name, payload, adapter_name, adapter_config, tcfg, opts, str(out_dir)) for name, payload in items]
    results: List[Dict[str, Any]] = []
    t0 = time.time()
    if opts.workers and opts.workers > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=opts.workers) as pool:
            futures = [pool.submit(_worker, j) for j in jobs]
            for i, fut in enumerate(as_completed(futures), 1):
                results.append(fut.result())
                if i % 25 == 0 or i == len(jobs):
                    log(f"  {i}/{len(jobs)} cells done ({time.time() - t0:.0f} s)")
    else:
        for i, j in enumerate(jobs, 1):
            results.append(_worker(j))
            if i % 25 == 0 or i == len(jobs):
                log(f"  {i}/{len(jobs)} cells done ({time.time() - t0:.0f} s)")

    results.sort(key=lambda r: r["path"])
    paths = [Path(r["path"]) for r in results]
    diag = pd.DataFrame([r["diag"] for r in results])
    diag.to_csv(out_dir / "tracking_diagnostics.csv", index=False)
    meta = {
        "snowtrack_version": __version__,
        "adapter": adapter_name,
        "adapter_config": adapter_config,
        "tracker_config": tcfg.__dict__,
        "options": {k: (list(v) if isinstance(v, tuple) else v) for k, v in opts.__dict__.items()},
        "source": str(source),
        "n_cells": len(paths),
    }
    (out_dir / "tracking_run.json").write_text(json.dumps(meta, indent=2, default=str), encoding="utf-8")

    if opts.summarize:
        frames = [pd.read_csv(p, parse_dates=["DATE"], compression="infer").assign(GRID_ID=_grid_id(p)) for p in paths]
        write_summaries(aggregate(frames), out_dir)
        log(f"  summaries written to {out_dir}")
    return paths


def _grid_id(p: Path) -> str:
    name = p.name
    for suf in (".gz", ".csv"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    return name[len("tracked_"):] if name.startswith("tracked_") else name
