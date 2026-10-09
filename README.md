# snowtrack

[![tests](https://github.com/rogerzzl/snow_track/actions/workflows/ci.yml/badge.svg?branch=v2.0)](https://github.com/rogerzzl/snow_track/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Where does the snowmelt go?** `snowtrack` takes the storages and fluxes that
any hydrological or land-surface model already writes and splits the runoff
of every cell, every time step, by *source* (rain, seasonal snowmelt, glacier
melt, ... any set) **and** by *pathway*: water that leaves through direct
surface runoff versus water that infiltrated, was stored (as liquid or soil
ice) and later left through drainage/baseflow. It is the dual-phase,
multi-layer source-tracking scheme of

> Zhu, Z., Tang, X., Zhang, J., Jin, J., Gao, C., Zhang, S., Xu, H., Tang, Y., Liu, Y., Guan, T.
> *Subsurface-routed snowmelt buffers river response to warming on the Tibetan Plateau.*

made model-agnostic (any number of soil layers, with or without an ice
phase, any set of sources) and packaged so that the results of the paper can
be reproduced from the raw model output with a few commands.

```
model output ──adapter──► canonical table ──tracker──► tracked_<cell>.csv ──route──► outlet series ──metrics──► T50, buffering, lag
 (VIC, CLM, Noah-MP,       S_RAIN S_SNOW ...            RUNOFF_SNOW_SURFACE        (Lohmann)                   Figure 1, Table S4,
  HBV, SWAT, your own)     W1_LIQ Q1_BASE ...           RUNOFF_SNOW_SUBSURFACE ...                            Figure 6 of the paper
```

The tracker does not modify the model: every bulk flux is partitioned among
sources in proportion to the composition of the store it leaves, and the
tracked storages are reconciled to the model's bulk storages every step, so
the decomposition is exactly mass-consistent ([`docs/method.md`](docs/method.md)).

## Quick start (10 minutes)

```bash
git clone https://github.com/rogerzzl/snow_track.git
cd snow_track
git checkout v2.0                     # version 2 lives on this branch (tag v2.0.0)
python -m pip install -e ".[dev]"      # Python >= 3.9; numpy, pandas, pyyaml (+ scipy, matplotlib, xarray for extras)
python -m pytest                       # ~1 min: unit tests + regression against the archived study outputs
snowtrack demo                         # synthetic cell, a real VIC-glacier cell, a basin outlet of the study
```

`snowtrack demo` ends with lines like

```
[1] synthetic 3-layer dual-phase cells: 2 tracked; max |closure error| = 9.9e-09 mm; snow share of runoff = 0.221 (...)
[2] VIC-glacier cell fluxes_31.0500_97.1500 (21915 days): max |diff| vs archived study output: surface 6.7e-16 mm, subsurface 1.0e-16 mm
[3] Changdu outlet 1968-2018: T50 advance (days/decade) surface snowmelt 8.73, subsurface snowmelt 1.53 (Table S4: 8.73 / 1.53); timing buffering index B = 7.95
```

Alternatives: `conda env create -f environment.yml`, or
`docker build -t snowtrack . && docker run --rm snowtrack` (runs the demo).
On Windows, `run_demo.ps1` does install + tests + demo.

## Use it with your model

1. Export, per cell and time step, what reaches the ground from each source,
   direct surface runoff, and for each store you want resolved its
   end-of-step storage and the fluxes leaving it (`snowtrack schema`,
   [`docs/input_schema.md`](docs/input_schema.md)). One store is enough;
   several layers and an ice phase are used when they are exported.
2. Tell `snowtrack` how to read it – pick an adapter
   ([`docs/adapters.md`](docs/adapters.md)):

   | your output | adapter | start from |
   |---|---|---|
   | VIC 4.x classic `fluxes_*` files | `vic-classic` | `configs/vic_classic_zhu2026.yaml` |
   | any CSV / text table per cell (HBV, SWAT, GR4J, ...) | `mapping` | `configs/mapping_template.yaml`, [`examples/04_bring_your_own_model`](examples/04_bring_your_own_model) |
   | gridded / HRU netCDF (CLM, Noah-MP, VIC image, SUMMA, ...) | `netcdf` | `configs/examples/clm_netcdf.yaml` |
   | already canonical | `canonical` | – |

3. Validate, run, route, analyse:

```bash
snowtrack validate -i model_out -c my_mapping.yaml                   # mapping check + water-budget residuals
snowtrack run      -i model_out -o tracked -c my_mapping.yaml --workers 8
snowtrack route    -i tracked -o outlet.csv --outlet 31.05 97.15 --surface-velocity <C_s> --base-velocity <C_b>
snowtrack metrics  --outlet outlet.csv -o metrics --start 1968 --end 2018
```

`run` writes `tracked_<cell>.csv` (snowmelt-derived surface and subsurface
runoff, the same for every other source, fractions, storage shares,
diagnostics) plus annual / monthly / overall summaries. `route` applies the
Lohmann channel routing scheme; the wave velocities of the surface and
subsurface pathways (and the diffusivity) are calibration parameters that you
calibrate against gauged discharge together with your model (`--cells` takes
flow-path distances or lags from your own routing instead). `metrics` writes
the monthly contribution climatology, annual T50 trends and buffering index,
late-season subsidy and apparent-lag diagnostics of the paper.

Python API:

```python
from snowtrack import CellData, SourceTracker, TrackerConfig
out = SourceTracker(TrackerConfig()).run(CellData.from_dataframe(df))   # df in the canonical schema
```

## Reproduce the paper

See [`docs/reproduce_paper.md`](docs/reproduce_paper.md) for the full chain
from raw VIC-glacier output to the paper's metrics, runtimes and the
data-availability table.

| level | what | command | result |
|---|---|---|---|
| 1 | scheme properties (closure, reduction, conventions) | `python -m pytest` | ~50 tests, Linux + Windows CI |
| 2 | per-cell tracking from raw VIC-glacier files | `python examples/02_vic_classic_changdu/run_example.py` | four real cells match the archived study output to ~1e-15 mm |
| 3 | basin-scale metrics from a routed outlet series of the study (Changdu) | `python examples/03_basin_outlet_metrics/basin_outlet_metrics.py` | Figure-1 climatology, Table S4 trends, Figure 6 lag diagnostics |
| 4 | full basin runs | `snowtrack run ... --composition legacy --skip-input-years 2`, then `snowtrack route` with calibrated velocities, on the complete VIC-glacier output (on request, ~1.4 GB per basin) | ~5 min per basin on 8 cores |

`--composition legacy` selects the numerical convention of the archived
study run, which the general tracker reproduces to machine precision
(`tests/test_regression_vic.py`).

## Repository layout

```
src/snowtrack/        package: schema.py (input), tracker.py (scheme), adapters/, pipeline.py, routing.py,
                      metrics.py, lagdiag.py, legacy.py (literal transcription of the archived code), cli.py
configs/              adapter configurations (VIC of the study, mapping template, CLM / Noah-MP / HBV examples)
examples/01..04       minimal dual-phase cell · real VIC-glacier cells + archived outputs · basin outlet metrics · your own model
docs/                 method.md (equations ↔ code), input_schema.md, adapters.md, reproduce_paper.md, faq.md
tests/                unit, regression and CLI tests (run by GitHub Actions on every push)
```

## Citing

Please cite the article above and the software (`CITATION.cff`; a Zenodo DOI
is minted for each tagged release). Contact: Zhanliang Zhu
(rogerzzl1997@gmail.com), Xiongpeng Tang (xptang@bnu.edu.cn).

## License

MIT (code and the bundled example data derived from the study's model runs).
