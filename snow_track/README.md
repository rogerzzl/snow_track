# Snow Track Checked

This folder contains a strict implementation of the dual-phase, three-layer
source-tracking scheme described in Section S3 of `supplementary_information.md`.

## Input columns

Each grid file must be CSV or whitespace-delimited and include:

- `YEAR`, `MONTH`, `DAY`
- `S_RAIN`, `S_SNOW`, `S_GLACIER`
- `INFIL`, `RUNOFF_SURF`
- For each layer `1..3`: `Wn_LIQ`, `Wn_ICE`, `Qn_PERC`, `Qn_BASE`,
  `En_TR`, `PHIn_FREEZE`, `PHIn_MELT`

Storages are in millimeters. Flux columns are interpreted as millimeters per
day and are multiplied by `--dt` during the explicit update. The sample file in
`sampledata/input/grid_001.csv` uses daily time steps with `--dt 1`.

## Run sample

```powershell
python .\snow_track.py --input-dir .\sampledata\input --output-dir .\sampledata\output --strict-surface-closure
```

## Run tests

```powershell
python .\test_snow_track.py
```
