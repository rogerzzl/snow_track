# Reproducing Zhu et al. – from raw model output to the paper's metrics

## The chain

| step | what | tool | output |
|---|---|---|---|
| 0 | VIC-glacier 4.2.d (classic driver, daily, 0.1°, 1961–2022, frozen-soil scheme), calibrated per basin against ERA5-Land SWE and gauged monthly discharge (Methods, S2) | VIC + glacier module, outside this repository | one `fluxes_<lat>_<lon>` file per cell with 14 columns |
| 1 | source tracking of every cell | `snowtrack run -c configs/vic_classic_zhu2026.yaml --composition legacy --skip-input-years 2` | `tracked_fluxes_*.csv` (1963–2022) |
| 2 | Lohmann routing of the tracked components to the gauge, with the calibrated wave velocities of the surface and subsurface pathways | `snowtrack route --outlet LAT LON --surface-velocity <C_s> --base-velocity <C_b>` | `<basin>_snowmelt.csv` |
| 3 | metrics | `snowtrack metrics`, `snowtrack.metrics`, `snowtrack.lagdiag` | fractions, T50, buffering index, late-season subsidy, apparent lag |

## What is in this repository

| data | status | where |
|---|---|---|
| raw VIC-glacier output, four cells of the Changdu basin (outlet, snow-richest, median, rain-dominated) | included, unmodified | `examples/02_vic_classic_changdu/input/` |
| archived per-cell tracking output for the same four cells | included | `examples/02_vic_classic_changdu/expected/` |
| routed snowmelt pathway series at the Changdu gauge | included | `examples/03_basin_outlet_metrics/data/` |
| complete raw VIC-glacier output (603 cells for Changdu, ~1.4 GB per basin), complete per-cell tracking output, routed series of the other basins | on request from the authors (too large for GitHub) | – |
| VIC-glacier parameter files, forcing (CN05.1), ERA5-Land SWE, gauge discharge | described in Methods / S1–S2 with their original sources; gauge discharge is not public | – |

## Commands and expected results

```bash
python -m pip install -e ".[dev]"
python -m pytest                                              # level 1: tests incl. regression against the archive
python examples/02_vic_classic_changdu/run_example.py         # level 2
python examples/03_basin_outlet_metrics/basin_outlet_metrics.py   # level 3
```

Level 2 prints, for each of the four cells, the maximum absolute difference
between the new output and the archived output of the study (surface and
subsurface snowmelt runoff): of the order of 1e-15 mm with
`--composition legacy`. The default convention (`end`) gives identical
surface runoff and annual subsurface snowmelt within 2 %.

Level 3 prints the Changdu metrics: monthly contribution climatology
(Figure 1), T50 advance of surface and subsurface snowmelt (8.73 and 1.53 days
per decade, Table S4), timing buffering index, late-season subsidy
(Figure 5c) and the apparent-lag diagnostics (Figure 6).

Runtime: the test suite ~1 min; level 2 ~15 s; level 3 ~10 s; a full basin
(603 cells, 62 years) ~3 s per cell and core, i.e. ~5 min with 8 workers.

## Settings of the archived run

* **Tracker**: adapter configuration `configs/vic_classic_zhu2026.yaml`
  (OUTVAR list of the study, `S_RAIN = OUT_RAINF`, `S_SNOW = OUT_SNOW_MELT`,
  `RUNOFF_SURF = OUT_RUNOFF`, baseflow `OUT_BASEFLOW`, soil evapotranspiration
  `OUT_EVAP` minus sublimation from the snow-water-equivalent balance);
  `--composition legacy` (withdrawals carry the start-of-step composition,
  the reported baseflow the end-of-step composition, storages renormalised
  to the bulk storage every day); the first two input years (1961–1962)
  dropped before tracking (`--skip-input-years 2`), after which the share of
  antecedent (UNKNOWN) water decays to < 0.1 % by 1968 (Figure S6).
* **Routing**: Lohmann channel routing from every cell to the gauge
  (`docs/method.md`, section 9). The wave velocities of the surface and
  subsurface pathways and the diffusivity are calibration parameters of the
  study; pass your own calibrated values to `snowtrack route`.
* **Table S4**: calendar-year T50 of surface snowmelt, subsurface snowmelt
  and total runoff at the gauge, OLS slope 1968–2018, days per decade,
  advance positive; B = trend_surf − trend_total (`snowtrack.metrics.buffering_index`).
* **S8**: July–October window, seasonal totals below 1 mm excluded, minimum
  30-day moving-average runoff, early 1968–1988 vs late 1998–2018
  (`late_season_lowflow`).
* **S9**: melt-response years April to March; T50 only for years whose
  volume exceeds the 5th percentile of the annual volumes; lag bins 0–30,
  30–90, 90–180, 180–365, 365–730, 730–1825 days; non-negative least squares
  (`snowtrack.lagdiag`).

`snowtrack.legacy` is a literal transcription of the archived tracking code;
`tests/test_regression_vic.py` shows that it, the general tracker in
`legacy` mode and the archived files agree to machine precision.

## Reproducing from scratch (full basins)

1. Obtain the VIC-glacier output of a basin (contact the authors), or run
   VIC-glacier yourself with the configuration of S2 and the OUTVAR list of
   `configs/vic_classic_zhu2026.yaml`.
2. `snowtrack run -i <vic_out> -o <basin>_tracked -c configs/vic_classic_zhu2026.yaml --composition legacy --skip-input-years 2 --workers 8 --columns compact`
3. `snowtrack route -i <basin>_tracked -o <basin>_snowmelt.csv --outlet <gauge lat> <gauge lon> --surface-velocity <C_s> --base-velocity <C_b>`
   with the velocities calibrated against the gauged discharge (gauge cells:
   Changdu 31.05 97.15, Dajin 31.55 102.05, Daofu 31.05 101.05, Ganzi 31.65
   99.95, Jiayuqiao 30.85 96.25, Kaqun 37.95 76.85, Lazi 29.15 87.55,
   Tangnaihai 35.45 100.05, Xiangtang 36.35 102.85, Zhimenda 33.05 97.25).
4. `snowtrack metrics --outlet <basin>_snowmelt.csv -o <basin>_metrics --start 1968 --end 2018`.
