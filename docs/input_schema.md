# Canonical input schema

One table per model cell (grid cell, HRU, sub-basin or lumped catchment), one
row per time step, sorted in time. Column names are case-insensitive. Print
this table with `snowtrack schema`.

| column | meaning | status |
|---|---|---|
| `DATE` **or** `YEAR, MONTH, DAY` (`HOUR`) | time stamp of the step | required |
| `S_RAIN`, `S_SNOW`, `S_GLACIER`, `S_<ANY>` | water reaching the ground surface from each source during the step | at least one |
| `RUNOFF_SURF` | direct surface runoff leaving the cell during the step | required |
| `INFIL` | infiltration into layer 1 during the step | optional, default `max(sum(S_*) - RUNOFF_SURF, 0)` |
| `W<n>_LIQ` | liquid storage of layer n at the **end** of the step | `W1_LIQ` required; layers must be numbered 1..N |
| `W<n>_ICE` | ice storage of layer n at the end of the step | optional (0) |
| `Q<n>_PERC` | percolation from layer n into layer n+1 (loss from the column for n = N) | optional (0) |
| `Q<n>_BASE` | drainage / baseflow from layer n to the channel | optional (0); at least one layer should have it |
| `E<n>_TR` | evapotranspiration withdrawn from layer n | optional (0) |
| `PHI<n>_FREEZE` | liquid-to-ice phase change in layer n | optional; inferred from `W<n>_ICE` when absent |
| `PHI<n>_MELT` | ice-to-liquid phase change in layer n | optional; inferred from `W<n>_ICE` when absent |

Units: storages in mm (or any depth unit), fluxes as depth **per time step**
in the same unit (mm per day for daily output). Sub-daily output works with a
`DATE` column or an `HOUR` column. Other columns are ignored (pass them to the
output with `--keep COL`).

## Semantics that matter

* **Sources are what reaches the ground.** `S_SNOW` is melt leaving the
  snowpack (plus rain-on-snow drainage if the model routes rain through the
  pack), not snowfall. `S_RAIN` is rain reaching the ground (after canopy
  interception if the model has a canopy). The tracker never uses snowfall or
  snow water equivalent.
* **`RUNOFF_SURF` is the surface pathway**, everything that leaves the cell
  without entering the stores you list. If your model has no infiltration
  term, leave `INFIL` out.
* **Stores are what you want resolved.** One store (`W1_LIQ` = total soil +
  groundwater storage, `Q1_BASE` = total baseflow) is the minimum. Several
  stores need the fluxes between them (`Q<n>_PERC`) to be exported too.
* **Evapotranspiration**: only the part taken from the listed stores.
  Canopy evaporation and snow sublimation do not come from the soil and must
  not be included (the VIC adapter subtracts sublimation from the bulk
  evapotranspiration for this reason).
* **End-of-step storages.** If your model writes storages at the start of the
  step, shift them by one step.

## Water-budget residuals

`snowtrack validate` reports, per cell, the mean absolute residual of

```
layer n:   [W_n(t) - W_n(t-1)] - [inflow_n - Q_n^perc - Q_n^base - E_n^tr]        (W = liquid + ice)
surface:   sum(S_*) - INFIL - RUNOFF_SURF
```

Zero means the exported fluxes close the exported storages; the tracker then
needs no reconciliation. Non-zero residuals are the fluxes your model did not
export (e.g. lateral flow, canopy drip, numerical closure). The tracker
absorbs them by rescaling the tracked storages to the bulk storage every
step and reports them per step (`BUDGET_RESIDUAL_LIQ<n>`). Large residuals
relative to the fluxes are a sign of a wrong mapping – check units and the
start/end-of-step convention first.

## Output columns (`tracked_<cell>.csv`)

| column | meaning |
|---|---|
| `TOTAL_RUNOFF`, `RUNOFF_SURF`, `BASEFLOW` | bulk runoff of the step (`RUNOFF_SURF + sum Q<n>_BASE`) |
| `RUNOFF_<SRC>_SURFACE` | part of `RUNOFF_SURF` originating from source SRC |
| `RUNOFF_<SRC>_SUBSURFACE` | part of `BASEFLOW` originating from source SRC |
| `RUNOFF_<SRC>_TOTAL` | their sum |
| `F_<SRC>_SURFACE / _SUBSURFACE / _TOTAL` | the same as fractions of `TOTAL_RUNOFF` (NaN when there is no runoff) |
| `FW_<SRC>` | share of SRC in the total liquid storage of the column |
| `W<n>_<SRC>_LIQ`, `W<n>_<SRC>_ICE` | source-resolved storages (`--columns full`) |
| `BUDGET_RESIDUAL_LIQ<n>`, `_ICE<n>` | bulk storage minus tracked storage before reconciliation |
| `CLIPPED_MASS`, `UNKNOWN_INJECTED` | mass clipped at zero / assigned to UNKNOWN during the step |
| `RUNOFF_SOURCE_CLOSURE_ERROR`, `W<n>_*_CLOSURE_ERROR` | closure checks (zero to round-off) |

`SRC` runs over the sources found in the input plus `UNKNOWN`.
`--columns compact` keeps only the time and runoff component columns.
