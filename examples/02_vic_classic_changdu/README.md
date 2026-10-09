# Example 02 – real VIC-glacier grid cells of the study (Changdu basin)

`input/` holds four unmodified daily output files of the calibrated VIC-glacier
run used in the paper (classic driver, 1961–2022, 14 columns, see
`vic_classic_zhu2026.yaml`):

| file | cell | why it is here |
|---|---|---|
| `fluxes_31.0500_97.1500.gz` | basin outlet cell (Changdu gauge) | outlet |
| `fluxes_32.1500_94.9500.gz` | highest snowmelt share of runoff in the basin (~14 %) | snow-rich headwater |
| `fluxes_31.7500_96.1500.gz` | median snowmelt share (~5 %) | typical cell |
| `fluxes_32.0500_96.4500.gz` | lowest snowmelt share (~1 %) | rain-dominated cell |

`expected/` holds, for the same cells, the per-cell tracking output that was
produced for the paper (`tracked_fluxes_*`, 1963–2022, archived convention).

## Reproduce

```bash
python examples/02_vic_classic_changdu/run_example.py
```

or by hand:

```bash
snowtrack run -i examples/02_vic_classic_changdu/input -o examples/02_vic_classic_changdu/output \
              -c examples/02_vic_classic_changdu/vic_classic_zhu2026.yaml \
              --composition legacy --skip-input-years 2 --workers 4
```

The script compares the new output with `expected/` (maximum absolute
difference of the order of 1e-15 mm) and then runs the default scheme
(`--composition end`) to show how much the choice of mixing convention
matters (annual subsurface snowmelt differs by < 1 %).

To run the full basin, point `-i` to the directory with all 603 `fluxes_*`
files of the VIC-glacier run (available from the authors on request; ~1.4 GB)
and route the result to the gauge with the Lohmann scheme:

```bash
snowtrack run   -i <vic_out> -o changdu_tracked -c vic_classic_zhu2026.yaml --composition legacy --skip-input-years 2 --workers 8
snowtrack route -i changdu_tracked -o changdu_outlet.csv --outlet 31.05 97.15 --surface-velocity <C_s> --base-velocity <C_b>
snowtrack metrics --outlet changdu_outlet.csv -o changdu_metrics --start 1968 --end 2018
```

The wave velocities of the surface and subsurface pathways are calibration
parameters (calibrated against the gauged discharge together with the
hydrological model). Runtime: about 3 s per cell and core.
