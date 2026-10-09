# Example 03 – metrics of a routed basin-outlet series

`data/changdu_snowmelt.csv.gz` is the daily series at the Changdu gauge
(Lancang River, 1961–2022) obtained by tracking every VIC-glacier grid cell of
the basin (Example 02) and routing the components to the gauge: snowmelt-
derived surface runoff, snowmelt-derived subsurface runoff and total runoff.
It is one of the ten basin series behind Figures 1–6 and Table S4 of the
paper and serves as the input example for `snowtrack metrics`. Columns follow
the archived naming (`RUNOFF_SNOW_SURFACE_OUTLET`, `RUNOFF_SNOW_BASE_OUTLET`,
`TOTAL_RUNOFF_OUTLET`); the metrics are ratios and timings and do not depend
on the depth scaling of the series.

```bash
python examples/03_basin_outlet_metrics/basin_outlet_metrics.py
# or
snowtrack metrics --outlet examples/03_basin_outlet_metrics/data/changdu_snowmelt.csv.gz -o changdu_metrics --start 1968 --end 2018
```

The script writes to `output/`:

| file | content | paper |
|---|---|---|
| `monthly_fraction_climatology.csv` | monthly mean contribution (%) of surface and subsurface snowmelt to total runoff | Figure 1 |
| `annual_fractions.csv` | annual snowmelt pathway fractions | Figure 2 |
| `timing_trends.csv` | OLS trends of annual T50 (days per decade, advance positive) | Figure 3, Table S4 |
| `summary.json` | timing buffering index `B = trend_surf − trend_total`, late-season subsidy, apparent-lag summary | S7, S8, S9 |
| `late_season_per_year.csv` | July–October subsurface snowmelt fraction and 30-day low-flow level per year | Figure 5c |
| `distributed_lag_kernel.csv`, `annual_subsurface_delay.csv`, `annual_carryover.csv`, `next_winter_release.csv` | non-negative distributed-lag model, T50 delay of subsurface vs surface snowmelt per melt-response year, apparent origin-year allocation | Figure 6, S9 |
| `changdu_summary.png` | two-panel summary figure | – |

For Changdu the script prints the T50 advance of surface (8.73 days per
decade) and subsurface (1.53) snowmelt as listed in Table S4. Series for the
other basins and the complete per-cell outputs are available from the authors
on request; `snowtrack route` produces the same kind of series from any
tracked run (`docs/reproduce_paper.md`).
