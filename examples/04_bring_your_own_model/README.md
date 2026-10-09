# Example 04 – bring your own model

Connect any model in three steps. The toy model here is an HBV-like lumped
bucket model with its own column names; replace it with your model output.

```bash
cd examples/04_bring_your_own_model
python make_toy_model_output.py                 # 1. your model writes per-cell tables (here: model_output/*.csv)
snowtrack validate -i model_output -c mapping.yaml   # 2. check the mapping; see the water-budget residuals
snowtrack run -i model_output -o tracked -c mapping.yaml   # 3. track
```

`mapping.yaml` is the only thing you write. Each canonical column is a raw
column name or a pandas expression over raw columns (`"suz + slz"`). What you
need from your model:

| canonical | what it is | the toy model's name |
|---|---|---|
| `S_RAIN`, `S_SNOW` (`S_GLACIER`, …) | water reaching the ground from each source | `rain`, `snowmelt` |
| `RUNOFF_SURF` | direct surface runoff | `qdirect` |
| `W<n>_LIQ` | end-of-step liquid storage of each store you want to resolve | `sm`, `suz + slz` |
| `Q<n>_PERC` | flow from store n to store n+1 | `perc` |
| `Q<n>_BASE` | flow from store n to the channel (subsurface pathway) | `q_upper + q_lower` |
| `E<n>_TR` | evapotranspiration taken from store n | `aet` |
| `W<n>_ICE`, `PHI<n>_*` | soil ice and freeze/thaw, only if your model has them | – |

One store is enough (`W1_LIQ` = total soil/groundwater storage, `Q1_BASE` =
total baseflow); more layers and an ice phase are used when your model exports
them. Gridded netCDF output uses the `netcdf` adapter instead
(`configs/examples/clm_netcdf.yaml`). See `docs/adapters.md`.

Interpretation of the output (`tracked/tracked_catchment_A.csv`):
`RUNOFF_SNOW_SURFACE` is snowmelt leaving through the surface pathway,
`RUNOFF_SNOW_SUBSURFACE` snowmelt that infiltrated and later left through the
response boxes, `F_SNOW_*` their shares of total runoff. `annual_source_fractions.csv`
and `monthly_source_fractions.csv` summarise per year / month.
