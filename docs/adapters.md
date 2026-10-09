# Connecting a hydrological model

`snowtrack` never touches a model's physics: it needs the model's exported
storages and fluxes in the canonical schema (`docs/input_schema.md`). An
*adapter* produces that schema. Four are built in; a few lines of Python add
another.

```
snowtrack validate -i <model output> -a <adapter> -c <config.yaml>   # check mapping + budget residuals
snowtrack convert  -i <model output> -o canonical -c <config.yaml>   # look at the canonical tables
snowtrack run      -i <model output> -o tracked   -c <config.yaml>   # track
```

| adapter | input | configuration |
|---|---|---|
| `canonical` | CSV / whitespace text already in the schema | none (`--adapter canonical`) |
| `vic-classic` | VIC 4.x classic-driver ASCII `fluxes_*` files | `configs/vic_classic_zhu2026.yaml` |
| `mapping` | any per-cell table (CSV/TSV, header or not) | `configs/mapping_template.yaml` |
| `netcdf` | gridded / HRU netCDF (CLM, Noah-MP, VIC image driver, SUMMA, ...) | `configs/examples/clm_netcdf.yaml` |

## VIC classic driver (`vic-classic`)

The classic driver writes one header-less file per cell whose columns follow
the `OUTVAR` list of the global parameter file after the time columns; list
them in the config (`outvars`, `nlayer`, `time_columns`). Per-layer variables
(`OUT_SOIL_MOIST`, `OUT_SOIL_LIQ`, `OUT_SOIL_ICE`, ...) expand into `nlayer`
columns. The adapter builds the canonical table

| canonical | VIC |
|---|---|
| `S_RAIN` | `OUT_RAINF` |
| `S_SNOW` | `OUT_SNOW_MELT` |
| `S_GLACIER` | only if your VIC variant writes glacier melt as an OUTVAR (`sources:` in the config) |
| `RUNOFF_SURF` | `OUT_RUNOFF` |
| `W1_LIQ` | sum of `OUT_SOIL_LIQ` layers, else sum of `OUT_SOIL_MOIST` (minus `OUT_SOIL_ICE` if written) |
| `W1_ICE` | sum of `OUT_SOIL_ICE` layers, if written (ice phase) |
| `Q1_BASE` | `OUT_BASEFLOW` |
| `E1_TR` | `OUT_EVAP` minus sublimation; sublimation from the SWE balance `SWE(t-1) + SNOWF - MELT - SWE(t)` or from `OUT_SUB_SNOW` if written |

Infiltration is `max(S_tot - RUNOFF, 0)` (tracker default). Soil ice is used
as a second phase when `OUT_SOIL_ICE` is among the OUTVARs. The VIC image
driver writes netCDF with the same variable names (`OUT_RAINF`, ...) – use the
`netcdf` adapter (`sum_over: nlayer`, or `isel` to keep layers apart).

## Any tabular output (`mapping`)

Write a YAML with the raw column (or a pandas expression) for each canonical
column; see `configs/mapping_template.yaml`, `configs/examples/lumped_hbv_csv.yaml`
and Example 04 for a two-store conceptual model (soil box + response boxes).
Typical mappings:

| model family | `S_RAIN` | `S_SNOW` | `RUNOFF_SURF` | stores |
|---|---|---|---|---|
| HBV-type | rain | snowmelt | 0 or infiltration-excess term | `W1_LIQ` = SM, `Q1_PERC` = recharge; `W2_LIQ` = SUZ + SLZ, `Q2_BASE` = Q0 + Q1 + Q2 |
| SWAT (HRU output) | PRECIP − SNOFALL | SNOMELT | SURQ | `W1_LIQ` = SW (soil), `Q1_PERC` = PERC, `E1_TR` = ET; `W2_LIQ` = shallow aquifer, `Q2_BASE` = GW_Q (+ LATQ as `Q1_BASE`) |
| GR4J / GR-family | rain | CemaNeige melt | Pn − Ps excess | `W1_LIQ` = production store, `Q1_PERC` = perc; `W2_LIQ` = routing store, `Q2_BASE` = Qr + Qd |
| Lumped bucket with one store | rain | melt | quick flow | `W1_LIQ` = store, `Q1_BASE` = slow flow |

Rules of thumb: a flux that leaves a store *to the channel* is `Q<n>_BASE`
(subsurface pathway); a flux from store n to store n+1 is `Q<n>_PERC`; a flux
to the channel that never entered a listed store is `RUNOFF_SURF`. Accumulated
variables must be differenced first; per-second units need `scale: 86400`.

## Gridded netCDF (`netcdf`)

```yaml
adapter: netcdf
time_dim: time
cell_dims: [lat, lon]          # or [hru], [gridcell], [south_north, west_east]
variables:
  S_RAIN: {var: RAIN, scale: 86400}
  S_SNOW: {var: QSNOMELT, scale: 86400}
  RUNOFF_SURF: {var: QOVER, scale: 86400}
  W1_LIQ: {var: SOILLIQ, sum_over: levsoi}
  W1_ICE: {var: SOILICE, sum_over: levsoi}
  Q1_BASE: {var: QDRAI, scale: 86400}
  E1_TR:  {var: QVEGT, scale: 86400}
```

Every cell of the grid becomes one tracked file (`cell_<lat>_<lon>`); cells
whose `mask_var` is all-NaN are skipped. A layered configuration selects
level indices with `isel` (e.g. `W1_LIQ: {var: SOILLIQ, isel: {levsoi: [0, 1, 2]}, sum_over: levsoi}`)
– but only if the model also exports the fluxes between those groups of
levels. Reading uses `xarray` (`pip install "snowtrack[netcdf]"`); multi-file
datasets are opened with `open_mfdataset` when a directory is given.

## Your own adapter

```python
from snowtrack.adapters import Adapter, register_adapter
import pandas as pd

@register_adapter
class MyModelAdapter(Adapter):
    name = "mymodel"
    patterns = ("*.out",)
    def load(self, path):
        raw = pd.read_fwf(path)
        return pd.DataFrame({"DATE": raw["date"], "S_RAIN": raw["rain"], "S_SNOW": raw["melt"],
                             "RUNOFF_SURF": raw["qs"], "W1_LIQ": raw["store"], "Q1_BASE": raw["qb"]})
```

Then `snowtrack run -a mymodel ...` (import your module first, e.g. via a small
driver script that calls `snowtrack.cli.main([...])`), or use the Python API
directly: `SourceTracker(TrackerConfig()).run(CellData.from_dataframe(df))`.

## Checklist before trusting a mapping

1. `snowtrack validate` shows the sources, layers and ice flag you expect.
2. The layer residuals are small compared with the fluxes (see
   `docs/input_schema.md`).
3. `RUNOFF_UNKNOWN_*` decays within a few years in `tracked_*.csv`; if it does
   not, water enters the stores through a path you have not mapped.
4. `CLIPPED_MASS` is ~0; otherwise withdrawals exceed the stores (wrong store,
   start-of-step storages, or unit mismatch).
