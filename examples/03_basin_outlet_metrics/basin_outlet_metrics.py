"""Metrics of the paper for the bundled Changdu outlet series (Figures 1-6, Table S4, S7-S9)."""

import json
from pathlib import Path

import pandas as pd

from snowtrack.metrics import (annual_fractions, buffering_index, late_season_lowflow, monthly_fraction_climatology,
                               read_outlet_series, standardise_outlet_columns, timing_trends)

HERE = Path(__file__).resolve().parent
BASIN = "changdu"
TABLE_S4 = {"trend_surf": 8.73, "trend_sub": 1.53}      # Changdu row, days per decade, advance positive


def main() -> None:
    out = HERE / "output"
    out.mkdir(exist_ok=True)
    std = standardise_outlet_columns(read_outlet_series(HERE / "data" / f"{BASIN}_snowmelt.csv.gz"))
    std = std[(std.index.year >= 1968) & (std.index.year <= 2018)]

    clim = monthly_fraction_climatology(std)
    clim.to_csv(out / "monthly_fraction_climatology.csv")
    annual_fractions(std).to_csv(out / "annual_fractions.csv")
    tt = timing_trends(std)
    tt.to_csv(out / "timing_trends.csv")
    bi = buffering_index(std)
    ls = late_season_lowflow(std)
    ls["per_year"].to_csv(out / "late_season_per_year.csv")
    summary = {"basin": BASIN, "buffering_index": bi, "late_season": {k: v for k, v in ls.items() if k != "per_year"}}

    fit = None
    try:
        from snowtrack.lagdiag import allocate_origin_year, annual_subsurface_delay, fit_distributed_lag, next_winter_release
        fit = fit_distributed_lag(std["Q_SURF"].to_numpy(), std["Q_SUB"].to_numpy())
        fit["kernel"].to_csv(out / "distributed_lag_kernel.csv", index=False)
        delay = annual_subsurface_delay(std, 4, 1968, 2017, threshold="p5")
        delay.to_csv(out / "annual_subsurface_delay.csv")
        allocate_origin_year(std.index, std["Q_SURF"].to_numpy(), std["Q_SUB"].to_numpy(), fit["daily_kernel"]).to_csv(out / "annual_carryover.csv")
        next_winter_release(std).to_csv(out / "next_winter_release.csv")
        summary["apparent_lag"] = fit["summary"]
        summary["median_T50_delay_days"] = float(delay["delay_T50_days"].median())
    except ImportError:
        summary["apparent_lag"] = "skipped (scipy not installed)"
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")

    pd.set_option("display.width", 160)
    print(f"{BASIN}: monthly mean contribution to total runoff (%)")
    print(clim[["F_SURF_MEAN_PCT", "F_SUB_MEAN_PCT"]].round(2).T.to_string())
    print(f"\nT50 advance 1968-2018 (days per decade): surface snowmelt {bi['trend_surf']:.2f} (Table S4: {TABLE_S4['trend_surf']}), "
          f"subsurface snowmelt {bi['trend_sub']:.2f} (Table S4: {TABLE_S4['trend_sub']}), total runoff {bi['trend_total']:.2f}; "
          f"timing buffering index B = {bi['B']:.2f}")
    print(f"late season (Jul-Oct): subsurface snowmelt fraction early {ls['F_sub_early_pct']:.2f} % -> late {ls['F_sub_late_pct']:.2f} %; "
          f"30-day low flow change {ls['dQ_low']:.3f}")
    if fit is not None:
        s = fit["summary"]
        print(f"apparent lag: median {s['apparent_median_lag_days']:.0f} d, tau90 {s['apparent_tau90_days']:.0f} d, "
              f"F(>1 yr) {s['F_gt365d']:.2f}, model R2 {s['model_R2']:.2f}; median T50 delay of subsurface vs surface "
              f"{summary['median_T50_delay_days']:.0f} d")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
        ax[0].bar(clim.index - 0.2, clim["F_SURF_MEAN_PCT"], 0.4, label="surface snowmelt")
        ax[0].bar(clim.index + 0.2, clim["F_SUB_MEAN_PCT"], 0.4, label="subsurface snowmelt")
        ax[0].set_xlabel("month"); ax[0].set_ylabel("% of total runoff"); ax[0].legend(frameon=False)
        ax[0].set_title(f"{BASIN}: monthly mean contribution")
        if fit is not None:
            ax[1].bar(fit["kernel"]["age_bin"], fit["kernel"]["apparent_fraction"])
            ax[1].set_ylabel("apparent lag fraction"); ax[1].set_title("distributed-lag kernel (S9)")
            ax[1].tick_params(axis="x", rotation=45)
        fig.tight_layout(); fig.savefig(out / f"{BASIN}_summary.png", dpi=150); plt.close(fig)
    except ImportError:
        pass
    print(f"\noutputs in {out}")


if __name__ == "__main__":
    main()
