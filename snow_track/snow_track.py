#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dual-phase, three-layer snowmelt source tracker.

This implementation follows Section S3 of supplementary_information.md:
source-resolved liquid and ice storages are tracked for three soil layers and
four sources: Rain, Snow, Glacier, and Unknown initial water.
"""

from __future__ import annotations

import argparse
import glob
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


SOURCES: Tuple[str, ...] = ("RAIN", "SNOW", "GLACIER", "UNKNOWN")
ACTIVE_SOURCES: Tuple[str, ...] = ("RAIN", "SNOW", "GLACIER")
LAYERS: Tuple[int, ...] = (1, 2, 3)


@dataclass(frozen=True)
class TrackerConfig:
    """Numerical controls for the explicit source-tracking update."""

    dt: float = 1.0
    epsilon: float = 1.0e-8
    strict_surface_closure: bool = False
    closure_tolerance: float = 1.0e-6


class DualPhaseSnowmeltTracker:
    """Run the S3 source-tracking scheme for VIC-glacier style grid files."""

    def __init__(
        self,
        input_dir: str | os.PathLike[str],
        output_dir: str | os.PathLike[str],
        config: TrackerConfig | None = None,
        spinup_years: int = 0,
    ) -> None:
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.config = config or TrackerConfig()
        self.spinup_years = spinup_years
        self.output_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def required_columns() -> List[str]:
        columns = [
            "YEAR",
            "MONTH",
            "DAY",
            "S_RAIN",
            "S_SNOW",
            "S_GLACIER",
            "INFIL",
            "RUNOFF_SURF",
        ]
        for layer in LAYERS:
            columns.extend(
                [
                    f"W{layer}_LIQ",
                    f"W{layer}_ICE",
                    f"Q{layer}_PERC",
                    f"Q{layer}_BASE",
                    f"E{layer}_TR",
                    f"PHI{layer}_FREEZE",
                    f"PHI{layer}_MELT",
                ]
            )
        return columns

    def load_data(self, file_path: str | os.PathLike[str]) -> pd.DataFrame:
        df = pd.read_csv(file_path, sep=r"\s+|,", engine="python")
        df.columns = [str(column).strip().upper() for column in df.columns]
        self._validate_input_columns(df, file_path)
        df["DATE"] = pd.to_datetime(df[["YEAR", "MONTH", "DAY"]])
        return df

    def process_grid_cell(self, file_path: str | os.PathLike[str]) -> pd.DataFrame:
        df = self.load_data(file_path)
        if df.empty:
            raise ValueError(f"{file_path} has no rows after spinup filtering.")

        n_steps = len(df)
        n_layers = len(LAYERS)
        n_sources = len(SOURCES)
        source_index = {source: idx for idx, source in enumerate(SOURCES)}

        liquid = np.zeros((n_layers, n_sources), dtype=float)
        ice = np.zeros((n_layers, n_sources), dtype=float)

        for layer_pos, layer in enumerate(LAYERS):
            liquid[layer_pos, source_index["UNKNOWN"]] = float(df.at[0, f"W{layer}_LIQ"])
            ice[layer_pos, source_index["UNKNOWN"]] = float(df.at[0, f"W{layer}_ICE"])

        records: List[Dict[str, float | int | pd.Timestamp]] = []

        for row_pos, row in df.iterrows():
            fluxes = self._zero_flux_record(row)
            if row_pos > 0:
                fluxes = self._advance_one_step(row, liquid, ice, source_index)

            self._reconcile_to_bulk(row, liquid, ice, phase="LIQ")
            self._reconcile_to_bulk(row, liquid, ice, phase="ICE")
            records.append(self._build_output_record(row, liquid, ice, fluxes))

        result = pd.DataFrame.from_records(records)
        return self._drop_spinup_outputs(result)

    def run(self, n_processes: int = 1, analyze: bool = True) -> None:
        file_paths = self._input_files()
        if not file_paths:
            raise FileNotFoundError(f"No input files found in {self.input_dir}.")

        if n_processes != 1:
            raise ValueError("Use n_processes=1 for deterministic validation output.")

        results = []
        for file_path in file_paths:
            result = self.process_grid_cell(file_path)
            output_path = self.output_dir / f"tracked_{Path(file_path).name}"
            result.to_csv(output_path, index=False)
            results.append(result.assign(GRID_ID=Path(file_path).name))

        if analyze:
            self.analyze_results(results)

    def analyze_results(self, results: Iterable[pd.DataFrame] | None = None) -> None:
        if results is None:
            result_files = sorted(self.output_dir.glob("tracked_*"))
            if not result_files:
                raise FileNotFoundError(f"No tracked result files found in {self.output_dir}.")
            frames = [
                pd.read_csv(path, parse_dates=["DATE"]).assign(
                    GRID_ID=path.name.replace("tracked_", "", 1)
                )
                for path in result_files
            ]
        else:
            frames = list(results)

        all_results = pd.concat(frames, ignore_index=True)
        annual = self._aggregate_fraction(all_results, ["GRID_ID", "YEAR"])
        monthly = self._aggregate_fraction(all_results, ["GRID_ID", "YEAR", "MONTH"])
        overall = self._overall_fraction(all_results)

        annual.to_csv(self.output_dir / "annual_snowmelt_fraction.csv", index=False)
        monthly.to_csv(self.output_dir / "monthly_snowmelt_fraction.csv", index=False)
        overall.to_csv(self.output_dir / "overall_source_fractions.csv", index=False)

    def _advance_one_step(
        self,
        row: pd.Series,
        liquid: np.ndarray,
        ice: np.ndarray,
        source_index: Dict[str, int],
    ) -> Dict[str, float]:
        self._validate_surface_closure(row)

        surface_fractions = self._surface_source_fractions(row)
        infiltration = {
            source: float(row["INFIL"]) * surface_fractions[source]
            for source in ACTIVE_SOURCES
        }
        surface_runoff = {
            source: float(row["RUNOFF_SURF"]) * surface_fractions[source] * self.config.dt
            for source in ACTIVE_SOURCES
        }

        baseflow = {source: 0.0 for source in SOURCES}
        percolation_from_above = np.zeros(len(SOURCES), dtype=float)

        for layer_pos, layer in enumerate(LAYERS):
            liquid_fraction = self._composition(liquid[layer_pos])

            if layer_pos == 0:
                inflow = np.zeros(len(SOURCES), dtype=float)
                for source in ACTIVE_SOURCES:
                    inflow[source_index[source]] = infiltration[source]
            else:
                inflow = percolation_from_above.copy()

            q_perc = float(row[f"Q{layer}_PERC"]) * liquid_fraction
            q_base = float(row[f"Q{layer}_BASE"]) * liquid_fraction
            evap = float(row[f"E{layer}_TR"]) * liquid_fraction

            liquid[layer_pos] = liquid[layer_pos] + (
                inflow - q_perc - q_base - evap
            ) * self.config.dt
            self._enforce_nonnegative(liquid[layer_pos])

            for source, idx in source_index.items():
                baseflow[source] += q_base[idx] * self.config.dt

            percolation_from_above = q_perc

        for layer_pos, layer in enumerate(LAYERS):
            freeze_total = float(row[f"PHI{layer}_FREEZE"])
            melt_total = float(row[f"PHI{layer}_MELT"])

            liquid_fraction = self._composition(liquid[layer_pos])
            ice_fraction = self._composition(ice[layer_pos])

            if liquid[layer_pos].sum() < self.config.epsilon:
                freeze_total = 0.0
            if ice[layer_pos].sum() < self.config.epsilon:
                melt_total = 0.0

            freeze = freeze_total * liquid_fraction * self.config.dt
            melt = melt_total * ice_fraction * self.config.dt

            liquid[layer_pos] = liquid[layer_pos] - freeze + melt
            ice[layer_pos] = ice[layer_pos] + freeze - melt

            self._enforce_nonnegative(liquid[layer_pos])
            self._enforce_nonnegative(ice[layer_pos])

        return {
            "RUNOFF_SURF_RAIN": surface_runoff["RAIN"],
            "RUNOFF_SURF_SNOW": surface_runoff["SNOW"],
            "RUNOFF_SURF_GLACIER": surface_runoff["GLACIER"],
            "BASEFLOW_RAIN": baseflow["RAIN"],
            "BASEFLOW_SNOW": baseflow["SNOW"],
            "BASEFLOW_GLACIER": baseflow["GLACIER"],
            "BASEFLOW_UNKNOWN": baseflow["UNKNOWN"],
        }

    def _surface_source_fractions(self, row: pd.Series) -> Dict[str, float]:
        supplies = {
            "RAIN": max(float(row["S_RAIN"]), 0.0),
            "SNOW": max(float(row["S_SNOW"]), 0.0),
            "GLACIER": max(float(row["S_GLACIER"]), 0.0),
        }
        total_supply = sum(supplies.values())
        if total_supply <= self.config.epsilon:
            return {source: 0.0 for source in ACTIVE_SOURCES}
        return {source: value / total_supply for source, value in supplies.items()}

    def _composition(self, storage_by_source: np.ndarray) -> np.ndarray:
        total = float(storage_by_source.sum())
        if total < self.config.epsilon:
            return np.zeros_like(storage_by_source)
        return storage_by_source / total

    def _reconcile_to_bulk(
        self,
        row: pd.Series,
        liquid: np.ndarray,
        ice: np.ndarray,
        phase: str,
    ) -> None:
        storage = liquid if phase == "LIQ" else ice
        for layer_pos, layer in enumerate(LAYERS):
            target = max(float(row[f"W{layer}_{phase}"]), 0.0)
            current = float(storage[layer_pos].sum())
            if target < self.config.epsilon:
                storage[layer_pos, :] = 0.0
            elif current < self.config.epsilon:
                storage[layer_pos, :] = 0.0
                storage[layer_pos, SOURCES.index("UNKNOWN")] = target
            else:
                storage[layer_pos, :] *= target / current

    def _enforce_nonnegative(self, values: np.ndarray) -> None:
        small_negative = (values < 0.0) & (values > -self.config.epsilon)
        values[small_negative] = 0.0
        if np.any(values < -self.config.epsilon):
            raise ValueError(
                "A source-resolved storage became negative. Check that input fluxes "
                "are consistent with bulk storages."
            )

    def _validate_surface_closure(self, row: pd.Series) -> None:
        supply = float(row["S_RAIN"] + row["S_SNOW"] + row["S_GLACIER"])
        routed = float(row["INFIL"] + row["RUNOFF_SURF"])
        residual = abs(supply - routed)
        if residual <= self.config.closure_tolerance:
            return
        message = (
            f"Surface closure residual on {row['DATE'].date()}: "
            f"S_tot={supply:.8f}, INFIL+RUNOFF_SURF={routed:.8f}."
        )
        if self.config.strict_surface_closure:
            raise ValueError(message)
        print(f"Warning: {message}")

    def _zero_flux_record(self, row: pd.Series) -> Dict[str, float]:
        return {
            "RUNOFF_SURF_RAIN": 0.0,
            "RUNOFF_SURF_SNOW": 0.0,
            "RUNOFF_SURF_GLACIER": 0.0,
            "BASEFLOW_RAIN": 0.0,
            "BASEFLOW_SNOW": 0.0,
            "BASEFLOW_GLACIER": 0.0,
            "BASEFLOW_UNKNOWN": 0.0,
        }

    def _build_output_record(
        self,
        row: pd.Series,
        liquid: np.ndarray,
        ice: np.ndarray,
        fluxes: Dict[str, float],
    ) -> Dict[str, float | int | pd.Timestamp]:
        surface_snow = fluxes["RUNOFF_SURF_SNOW"]
        subsurface_snow = fluxes["BASEFLOW_SNOW"]
        total_snow = surface_snow + subsurface_snow
        total_runoff = (float(row["RUNOFF_SURF"]) + sum(
            float(row[f"Q{layer}_BASE"]) for layer in LAYERS
        )) * self.config.dt

        record: Dict[str, float | int | pd.Timestamp] = {
            "DATE": row["DATE"],
            "YEAR": int(row["YEAR"]),
            "MONTH": int(row["MONTH"]),
            "DAY": int(row["DAY"]),
            "TOTAL_RUNOFF": total_runoff,
            "RUNOFF_SNOW_TOTAL": total_snow,
            "RUNOFF_SNOW_SURFACE": surface_snow,
            "RUNOFF_SNOW_SUBSURFACE": subsurface_snow,
            "F_SNOW_TOTAL": _safe_divide(total_snow, total_runoff),
            "F_SNOW_SURFACE": _safe_divide(surface_snow, total_runoff),
            "F_SNOW_SUBSURFACE": _safe_divide(subsurface_snow, total_runoff),
        }

        for source in SOURCES:
            source_total = (
                fluxes.get(f"RUNOFF_SURF_{source}", 0.0)
                + fluxes.get(f"BASEFLOW_{source}", 0.0)
            )
            record[f"RUNOFF_{source}_TOTAL"] = source_total
            record[f"RUNOFF_{source}_SURFACE"] = fluxes.get(f"RUNOFF_SURF_{source}", 0.0)
            record[f"RUNOFF_{source}_SUBSURFACE"] = fluxes.get(f"BASEFLOW_{source}", 0.0)

        for layer_pos, layer in enumerate(LAYERS):
            for source_pos, source in enumerate(SOURCES):
                record[f"W{layer}_{source}_LIQ"] = liquid[layer_pos, source_pos]
                record[f"W{layer}_{source}_ICE"] = ice[layer_pos, source_pos]
            record[f"W{layer}_LIQ_CLOSURE_ERROR"] = (
                sum(record[f"W{layer}_{source}_LIQ"] for source in SOURCES)
                - float(row[f"W{layer}_LIQ"])
            )
            record[f"W{layer}_ICE_CLOSURE_ERROR"] = (
                sum(record[f"W{layer}_{source}_ICE"] for source in SOURCES)
                - float(row[f"W{layer}_ICE"])
            )

        record["RUNOFF_SOURCE_CLOSURE_ERROR"] = (
            sum(record[f"RUNOFF_{source}_TOTAL"] for source in SOURCES) - total_runoff
        )
        return record

    def _drop_spinup_outputs(self, result: pd.DataFrame) -> pd.DataFrame:
        if self.spinup_years <= 0:
            return result
        first_year = int(result["YEAR"].iloc[0])
        first_analysis_date = pd.Timestamp(first_year + self.spinup_years, 1, 1)
        return result.loc[result["DATE"] >= first_analysis_date].reset_index(drop=True)

    def _validate_input_columns(
        self,
        df: pd.DataFrame,
        file_path: str | os.PathLike[str],
    ) -> None:
        missing = [column for column in self.required_columns() if column not in df.columns]
        if missing:
            missing_text = ", ".join(missing)
            raise ValueError(f"{file_path} is missing required columns: {missing_text}")

        numeric_columns = [column for column in self.required_columns() if column not in {"YEAR", "MONTH", "DAY"}]
        for column in numeric_columns:
            df[column] = pd.to_numeric(df[column], errors="raise")

    def _input_files(self) -> List[Path]:
        patterns = ["*.csv", "*.txt", "*.dat"]
        files: List[Path] = []
        for pattern in patterns:
            files.extend(Path(path) for path in glob.glob(str(self.input_dir / pattern)))
        return sorted(set(files))

    @staticmethod
    def _aggregate_fraction(df: pd.DataFrame, group_cols: List[str]) -> pd.DataFrame:
        rows = []
        for keys, group in df.groupby(group_cols, dropna=False):
            if not isinstance(keys, tuple):
                keys = (keys,)
            total_runoff = group["TOTAL_RUNOFF"].sum()
            row = dict(zip(group_cols, keys))
            row["RUNOFF_SNOW_TOTAL"] = group["RUNOFF_SNOW_TOTAL"].sum()
            row["TOTAL_RUNOFF"] = total_runoff
            row["F_SNOW_TOTAL"] = _safe_divide(row["RUNOFF_SNOW_TOTAL"], total_runoff)
            row["RUNOFF_SNOW_SURFACE"] = group["RUNOFF_SNOW_SURFACE"].sum()
            row["RUNOFF_SNOW_SUBSURFACE"] = group["RUNOFF_SNOW_SUBSURFACE"].sum()
            rows.append(row)
        return pd.DataFrame(rows)

    @staticmethod
    def _overall_fraction(df: pd.DataFrame) -> pd.DataFrame:
        total_runoff = df["TOTAL_RUNOFF"].sum()
        rows = []
        for source in SOURCES:
            source_total = df[f"RUNOFF_{source}_TOTAL"].sum()
            rows.append(
                {
                    "SOURCE": source,
                    "RUNOFF_SOURCE_TOTAL": source_total,
                    "TOTAL_RUNOFF": total_runoff,
                    "F_SOURCE_TOTAL": _safe_divide(source_total, total_runoff),
                }
            )
        return pd.DataFrame(rows)


def _safe_divide(numerator: float, denominator: float) -> float:
    if abs(float(denominator)) < 1.0e-12:
        return np.nan
    return float(numerator) / float(denominator)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the S3 dual-phase, three-layer snowmelt source tracker."
    )
    parser.add_argument("--input-dir", required=True, help="Directory with grid input files.")
    parser.add_argument("--output-dir", required=True, help="Directory for tracked outputs.")
    parser.add_argument("--spinup-years", type=int, default=0, help="Years to skip from the start.")
    parser.add_argument("--dt", type=float, default=1.0, help="Time step length in days.")
    parser.add_argument(
        "--strict-surface-closure",
        action="store_true",
        help="Fail if S_rain + S_snow + S_glacier differs from infiltration plus surface runoff.",
    )
    parser.add_argument(
        "--no-analyze",
        action="store_true",
        help="Skip annual, monthly, and overall aggregation files.",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    config = TrackerConfig(
        dt=args.dt,
        strict_surface_closure=args.strict_surface_closure,
    )
    tracker = DualPhaseSnowmeltTracker(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        config=config,
        spinup_years=args.spinup_years,
    )
    tracker.run(analyze=not args.no_analyze)


if __name__ == "__main__":
    main()
