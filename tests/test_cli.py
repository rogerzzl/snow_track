import json
import subprocess
import sys

import pandas as pd

from snowtrack.cli import main
from snowtrack.synthetic import make_synthetic_cell


def test_run_validate_route_metrics_roundtrip(tmp_path):
    inp = tmp_path / "in"; inp.mkdir()
    make_synthetic_cell(years=3, seed=1).to_csv(inp / "cell_30.0500_96.0500.csv", index=False)
    make_synthetic_cell(years=3, seed=2).to_csv(inp / "cell_30.1500_96.1500.csv", index=False)
    assert main(["validate", "-i", str(inp)]) == 0
    out = tmp_path / "tracked"
    assert main(["run", "-i", str(inp), "-o", str(out), "--workers", "2", "--columns", "standard"]) == 0
    files = sorted(out.glob("tracked_*.csv"))
    assert len(files) == 2
    assert (out / "annual_source_fractions.csv").exists() and (out / "domain_monthly_climatology.csv").exists()
    meta = json.loads((out / "tracking_run.json").read_text())
    assert meta["n_cells"] == 2
    outlet = tmp_path / "outlet.csv"
    assert main(["route", "-i", str(out), "-o", str(outlet), "--outlet", "30.05", "96.05",
                 "--surface-velocity", "1.5", "--base-velocity", "0.3"]) == 0
    df = pd.read_csv(outlet, comment="#")
    assert "RUNOFF_SNOW_SUBSURFACE_OUTLET" in df.columns and len(df) == 3 * 365
    mdir = tmp_path / "metrics"
    assert main(["metrics", "--outlet", str(outlet), "-o", str(mdir), "--early", "2001", "2001", "--late", "2003", "2003"]) == 0
    assert (mdir / "summary.json").exists() and (mdir / "timing_trends.csv").exists()


def test_convert_and_schema(tmp_path):
    inp = tmp_path / "in"; inp.mkdir()
    make_synthetic_cell(years=1, seed=1).to_csv(inp / "a.csv", index=False)
    assert main(["convert", "-i", str(inp), "-o", str(tmp_path / "canon")]) == 0
    assert (tmp_path / "canon" / "a.csv").exists()
    assert main(["schema", "--layers", "1"]) == 0


def test_demo_entry_point(tmp_path):
    proc = subprocess.run([sys.executable, "-m", "snowtrack.cli", "demo", "-o", str(tmp_path / "demo")],
                          capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stderr
    assert "[1] synthetic" in proc.stdout
