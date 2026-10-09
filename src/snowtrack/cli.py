"""Command-line interface: ``snowtrack <command> ...``."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from . import __version__
from .adapters import ADAPTERS, get_adapter, load_config
from .pipeline import RunOptions, run as run_pipeline
from .schema import CellData, SchemaError, canonical_column_summary
from .tracker import COMPOSITIONS, TrackerConfig

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = REPO_ROOT / "examples"


def _tracker_config(a) -> TrackerConfig:
    return TrackerConfig(
        composition=a.composition, strict=a.strict, strict_surface_closure=a.strict_surface_closure,
        spinup_years=a.spinup_years, store_storage_columns=(a.columns == "full"),
    )


def _run_options(a) -> RunOptions:
    return RunOptions(
        skip_input_years=a.skip_input_years, fill_missing=a.fill_missing, keep_columns=tuple(a.keep or ()),
        columns=a.columns, gzip=a.gzip, float_format=a.float_format, workers=a.workers,
        summarize=not a.no_summary, limit=a.limit,
    )


def _add_run_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--input-dir", "-i", required=True, help="directory of per-cell files (or a netCDF file)")
    p.add_argument("--output-dir", "-o", required=True)
    p.add_argument("--adapter", "-a", choices=sorted(ADAPTERS), default=None,
                   help="input format (default: 'adapter' key of --config, else canonical)")
    p.add_argument("--config", "-c", default=None, help="YAML adapter configuration (see configs/)")
    p.add_argument("--composition", choices=COMPOSITIONS, default="end",
                   help="composition carried by withdrawals: end (mix-then-withdraw, default), start, legacy (archived study run)")
    p.add_argument("--skip-input-years", type=int, default=0, help="drop the first N years of input before tracking")
    p.add_argument("--spinup-years", type=int, default=0, help="track from the start but drop the first N years from the output")
    p.add_argument("--workers", "-w", type=int, default=1, help="parallel processes")
    p.add_argument("--columns", choices=("full", "standard", "compact"), default="standard",
                   help="output columns: full (+source-resolved storages), standard, compact (runoff components only)")
    p.add_argument("--gzip", action="store_true", help="write tracked_*.csv.gz")
    p.add_argument("--float-format", default=None, help="e.g. %%.8g to shrink output files")
    p.add_argument("--fill-missing", type=float, default=None, help="replace NaN inputs by this value")
    p.add_argument("--keep", action="append", help="pass an input column through to the output (repeatable)")
    p.add_argument("--strict", action="store_true", help="fail when a tracked storage would go negative")
    p.add_argument("--strict-surface-closure", action="store_true", help="fail when sum(S_*) != INFIL + RUNOFF_SURF")
    p.add_argument("--no-summary", action="store_true", help="skip annual/monthly/overall summary tables")
    p.add_argument("--limit", type=int, default=None, help="only process the first N cells")


def cmd_run(a) -> int:
    cfg = load_config(a.config)
    paths = run_pipeline(a.input_dir, a.output_dir, a.adapter, cfg, _tracker_config(a), _run_options(a))
    print(f"wrote {len(paths)} tracked file(s) to {a.output_dir}")
    return 0


def cmd_convert(a) -> int:
    cfg = load_config(a.config)
    adapter = get_adapter(a.adapter, cfg)
    out = Path(a.output_dir); out.mkdir(parents=True, exist_ok=True)
    n = 0
    for name, payload in adapter.iter_items(a.input_dir):
        df = adapter.load(payload)
        df.to_csv(out / f"{name}.csv", index=False)
        n += 1
        if a.limit and n >= a.limit:
            break
    print(f"converted {n} cell(s) to canonical CSV in {out}")
    return 0


def cmd_validate(a) -> int:
    cfg = load_config(a.config)
    adapter = get_adapter(a.adapter, cfg)
    n = ok = 0
    for name, payload in adapter.iter_items(a.input_dir):
        n += 1
        try:
            cell = CellData.from_dataframe(adapter.load(payload), name=name, fill_missing=a.fill_missing)
            s = cell.summary()
            ok += 1
            print(f"[OK]   {name}: {s['steps']} steps {s['start']}..{s['end']}, sources={s['sources']}, "
                  f"layers={s['n_layers']}, ice={s['has_ice']}, mean|surface residual|={s['mean_abs_surface_residual']:.3g}, "
                  f"mean|layer residual|={['%.3g' % v for v in s['mean_abs_layer_residual']]}")
            if a.verbose:
                print("       provided:", json.dumps(s["provided"]))
        except (SchemaError, KeyError, ValueError) as exc:
            print(f"[FAIL] {name}: {exc}")
        if a.limit and n >= a.limit:
            break
    print(f"{ok}/{n} cell(s) valid")
    return 0 if ok == n else 1


def cmd_schema(a) -> int:
    pd.set_option("display.width", 160, "display.max_colwidth", 70)
    print(canonical_column_summary(a.layers).to_string(index=False))
    return 0


def cmd_route(a) -> int:
    from .routing import RoutingConfig, route, write_outlet
    cfg = RoutingConfig(surface_velocity=a.surface_velocity, base_velocity=a.base_velocity, diffusivity=a.diffusivity,
                        scheme=a.scheme, weight=a.weight, grid_resolution_deg=a.grid_resolution)
    table = pd.read_csv(a.cells) if a.cells else None
    outlet = tuple(a.outlet) if a.outlet else None
    df = route(a.tracked_dir, outlet, cfg, table, pattern=a.pattern)
    note = (f"scheme={a.scheme}; calibrated parameters: surface_velocity={a.surface_velocity} m/s, "
            f"base_velocity={a.base_velocity} m/s, diffusivity={a.diffusivity} m2/s")
    write_outlet(df, a.out, outlet, header_note=note)
    print(f"routed {df.attrs.get('n_cells')} cell(s) -> {a.out}  ({df.index.min().date()}..{df.index.max().date()})")
    return 0


def cmd_metrics(a) -> int:
    from .metrics import (annual_fractions, buffering_index, late_season_lowflow, monthly_fraction_climatology,
                          read_outlet_series, standardise_outlet_columns, timing_trends)
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    std = standardise_outlet_columns(read_outlet_series(a.outlet), a.source)
    if a.start:
        std = std[std.index.year >= a.start]
    if a.end:
        std = std[std.index.year <= a.end]
    monthly_fraction_climatology(std).to_csv(out / "monthly_fraction_climatology.csv")
    annual_fractions(std).to_csv(out / "annual_fractions.csv")
    tt = timing_trends(std, year_start_month=a.year_start_month)
    tt.to_csv(out / "timing_trends.csv")
    bi = buffering_index(std, year_start_month=a.year_start_month)
    ls = late_season_lowflow(std, early=tuple(a.early), late=tuple(a.late))
    ls["per_year"].to_csv(out / "late_season_per_year.csv")
    summary = {"buffering_index": bi, "late_season": {k: v for k, v in ls.items() if k != "per_year"}}
    if not a.no_lag:
        try:
            from .lagdiag import allocate_origin_year, annual_subsurface_delay, fit_distributed_lag, next_winter_release
            fit = fit_distributed_lag(std["Q_SURF"].to_numpy(), std["Q_SUB"].to_numpy())
            fit["kernel"].to_csv(out / "distributed_lag_kernel.csv", index=False)
            fit["daily_kernel"].to_csv(out / "distributed_lag_daily_kernel.csv", index=False)
            annual_subsurface_delay(std).to_csv(out / "annual_subsurface_delay.csv")
            allocate_origin_year(std.index, std["Q_SURF"].to_numpy(), std["Q_SUB"].to_numpy(), fit["daily_kernel"]).to_csv(out / "annual_carryover.csv")
            next_winter_release(std).to_csv(out / "next_winter_release.csv")
            summary["apparent_lag"] = fit["summary"]
        except ImportError as exc:
            summary["apparent_lag"] = f"skipped: {exc}"
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
    print(json.dumps(summary, indent=2, default=float))
    print(f"metrics written to {out}")
    return 0


def cmd_demo(a) -> int:
    """Synthetic closed cell + (if available) the bundled real examples."""
    from .synthetic import make_synthetic_cell
    out = Path(a.output_dir); out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print(f"snowtrack {__version__} demo -> {out}")

    # 1) synthetic three-layer dual-phase cell with an exactly closed budget
    syn_dir = out / "synthetic_input"; syn_dir.mkdir(exist_ok=True)
    make_synthetic_cell(years=4, seed=1).to_csv(syn_dir / "cell_A.csv", index=False)
    make_synthetic_cell(years=4, seed=2, glacier=False).to_csv(syn_dir / "cell_B.csv", index=False)
    paths = run_pipeline(syn_dir, out / "synthetic_tracked", "canonical", {}, TrackerConfig(), RunOptions(columns="full"), log=lambda *x: None)
    tr = pd.read_csv(paths[0])
    closure = [c for c in tr.columns if c.endswith("_CLOSURE_ERROR")]
    print(f"[1] synthetic 3-layer dual-phase cells: {len(paths)} tracked; max |closure error| = "
          f"{np.nanmax(np.abs(tr[closure].to_numpy())):.2e} mm; "
          f"snow share of runoff = {tr['RUNOFF_SNOW_TOTAL'].sum() / tr['TOTAL_RUNOFF'].sum():.3f} "
          f"(surface {tr['RUNOFF_SNOW_SURFACE'].sum() / tr['TOTAL_RUNOFF'].sum():.3f}, "
          f"subsurface {tr['RUNOFF_SNOW_SUBSURFACE'].sum() / tr['TOTAL_RUNOFF'].sum():.3f})")

    # 2) real VIC-glacier grid cell of the study, archived convention, regression check
    vic = EXAMPLES / "02_vic_classic_changdu"
    if vic.is_dir():
        cfg = load_config(vic / "vic_classic_zhu2026.yaml")
        paths = run_pipeline(vic / "input", out / "vic_tracked", "vic-classic", cfg,
                             TrackerConfig(composition="legacy", store_storage_columns=False),
                             RunOptions(skip_input_years=2, limit=1, summarize=False), log=lambda *x: None)
        got = pd.read_csv(paths[0], parse_dates=["DATE"])
        name = paths[0].name.replace("tracked_", "").replace(".csv", "")
        exp_path = vic / "expected" / f"tracked_{name}.csv.gz"
        if exp_path.exists():
            exp = pd.read_csv(exp_path, parse_dates=["DATE"])
            d_surf = np.nanmax(np.abs(got["RUNOFF_SNOW_SURFACE"].to_numpy() - exp["RUNOFF_SNOW_SURFACE"].to_numpy()))
            d_base = np.nanmax(np.abs(got["RUNOFF_SNOW_SUBSURFACE"].to_numpy() - exp["RUNOFF_SNOW_BASE"].to_numpy()))
            print(f"[2] VIC-glacier cell {name} ({len(got)} days): max |diff| vs archived study output: "
                  f"surface {d_surf:.1e} mm, subsurface {d_base:.1e} mm")
        else:
            print(f"[2] VIC-glacier cell {name} tracked ({len(got)} days)")
    else:
        print("[2] bundled VIC example not found (installed without the repository) - skipped")

    # 3) basin outlet metrics of the study
    f = EXAMPLES / "03_basin_outlet_metrics" / "data" / "changdu_snowmelt.csv.gz"
    if f.exists():
        from .metrics import buffering_index, read_outlet_series, standardise_outlet_columns
        std = standardise_outlet_columns(read_outlet_series(f))
        bi = buffering_index(std, 1968, 2018)
        print(f"[3] Changdu outlet 1968-2018: T50 advance (days/decade) surface snowmelt {bi['trend_surf']:.2f}, "
              f"subsurface snowmelt {bi['trend_sub']:.2f} (Table S4: 8.73 / 1.53); timing buffering index B = {bi['B']:.2f}")
    else:
        print("[3] bundled outlet example not found - skipped")
    print(f"demo finished in {time.time() - t0:.1f} s")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="snowtrack", description="Source tracking of snowmelt through surface and subsurface runoff pathways.")
    p.add_argument("--version", action="version", version=f"snowtrack {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("run", help="track every cell of a directory / dataset"); _add_run_args(s); s.set_defaults(func=cmd_run)

    s = sub.add_parser("convert", help="write canonical CSV files from any adapter")
    s.add_argument("--input-dir", "-i", required=True); s.add_argument("--output-dir", "-o", required=True)
    s.add_argument("--adapter", "-a", choices=sorted(ADAPTERS), default=None); s.add_argument("--config", "-c", default=None)
    s.add_argument("--limit", type=int, default=None); s.set_defaults(func=cmd_convert)

    s = sub.add_parser("validate", help="check inputs against the canonical schema and report budget residuals")
    s.add_argument("--input-dir", "-i", required=True); s.add_argument("--adapter", "-a", choices=sorted(ADAPTERS), default=None)
    s.add_argument("--config", "-c", default=None); s.add_argument("--fill-missing", type=float, default=None)
    s.add_argument("--limit", type=int, default=None); s.add_argument("--verbose", "-v", action="store_true"); s.set_defaults(func=cmd_validate)

    s = sub.add_parser("schema", help="print the canonical input columns"); s.add_argument("--layers", type=int, default=3); s.set_defaults(func=cmd_schema)

    s = sub.add_parser("route", help="route tracked cell components to a basin outlet (Lohmann channel routing)")
    s.add_argument("--tracked-dir", "-i", required=True); s.add_argument("--out", "-o", required=True)
    s.add_argument("--outlet", type=float, nargs=2, metavar=("LAT", "LON"), default=None, help="gauge location")
    s.add_argument("--surface-velocity", type=float, required=True, metavar="M_PER_S",
                   help="wave velocity of the surface pathway (calibration parameter)")
    s.add_argument("--base-velocity", type=float, required=True, metavar="M_PER_S",
                   help="wave velocity of the subsurface pathway (calibration parameter)")
    s.add_argument("--diffusivity", type=float, default=800.0, metavar="M2_PER_S",
                   help="diffusion-wave diffusivity of the Lohmann scheme (calibratable; default 800)")
    s.add_argument("--scheme", choices=("lohmann", "lag"), default="lohmann",
                   help="lohmann: diffusion-wave impulse response (default); lag: pure translation by distance / velocity")
    s.add_argument("--cells", default=None, help="CSV with cell_id and distance_m (or lag_surface, lag_base) [, weight|area]")
    s.add_argument("--weight", choices=("none", "area"), default="none"); s.add_argument("--grid-resolution", type=float, default=None)
    s.add_argument("--pattern", default="tracked_*"); s.set_defaults(func=cmd_route)

    s = sub.add_parser("metrics", help="timing, buffering, late-season and apparent-lag metrics of an outlet series")
    s.add_argument("--outlet", required=True, help="outlet CSV (from 'snowtrack route' or the bundled example)")
    s.add_argument("--out-dir", "-o", required=True); s.add_argument("--source", default="SNOW")
    s.add_argument("--start", type=int, default=None); s.add_argument("--end", type=int, default=None)
    s.add_argument("--year-start-month", type=int, default=1, help="1 = calendar year (T50 trends), 4 = melt-response year")
    s.add_argument("--early", type=int, nargs=2, default=(1968, 1988)); s.add_argument("--late", type=int, nargs=2, default=(1998, 2018))
    s.add_argument("--no-lag", action="store_true", help="skip the distributed-lag diagnostics (needs scipy)")
    s.set_defaults(func=cmd_metrics)

    s = sub.add_parser("demo", help="run the bundled examples end to end"); s.add_argument("--output-dir", "-o", default="snowtrack_demo"); s.set_defaults(func=cmd_demo)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
