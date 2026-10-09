"""Track the four bundled VIC-glacier cells and compare with the archived study output."""

from pathlib import Path

import numpy as np
import pandas as pd

from snowtrack.adapters import load_config
from snowtrack.pipeline import RunOptions, run
from snowtrack.tracker import TrackerConfig

HERE = Path(__file__).resolve().parent


def compare(out_dir: Path, label: str) -> None:
    print(f"\n{label}")
    print(f"{'cell':28s} {'max|dSurf| mm':>14s} {'max|dSub| mm':>13s} {'annual sub. rel. diff':>22s}")
    for p in sorted(out_dir.glob("tracked_*.csv")):
        name = p.name[len("tracked_"):-len(".csv")]
        got = pd.read_csv(p, parse_dates=["DATE"])
        exp = pd.read_csv(HERE / "expected" / f"tracked_{name}.csv.gz", parse_dates=["DATE"])
        d_s = np.nanmax(np.abs(got["RUNOFF_SNOW_SURFACE"].to_numpy() - exp["RUNOFF_SNOW_SURFACE"].to_numpy()))
        d_b = np.nanmax(np.abs(got["RUNOFF_SNOW_SUBSURFACE"].to_numpy() - exp["RUNOFF_SNOW_BASE"].to_numpy()))
        a_g = got.groupby("YEAR")["RUNOFF_SNOW_SUBSURFACE"].sum()
        a_e = exp.groupby("YEAR")["RUNOFF_SNOW_BASE"].sum()
        rel = float(np.nanmax(np.abs(a_g - a_e) / a_e.replace(0, np.nan)))
        print(f"{name:28s} {d_s:14.2e} {d_b:13.2e} {rel:22.2e}")


if __name__ == "__main__":
    cfg = load_config(HERE / "vic_classic_zhu2026.yaml")
    out_legacy = HERE / "output" / "legacy"
    run(HERE / "input", out_legacy, "vic-classic", cfg, TrackerConfig(composition="legacy"),
        RunOptions(skip_input_years=2, workers=4))
    compare(out_legacy, "archived convention (--composition legacy): should match to machine precision")

    out_end = HERE / "output" / "default"
    run(HERE / "input", out_end, "vic-classic", cfg, TrackerConfig(composition="end"),
        RunOptions(skip_input_years=2, workers=4))
    compare(out_end, "default convention (--composition end): same surface runoff, slightly different subsurface timing")

    ov = pd.read_csv(out_legacy / "overall_source_fractions.csv")
    print("\nOverall source shares of runoff, four cells pooled, 1963-2022:")
    print(ov[["SOURCE", "F_SOURCE_TOTAL", "F_SOURCE_SURFACE", "F_SOURCE_SUBSURFACE"]].round(4).to_string(index=False))
