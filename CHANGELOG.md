# Changelog

## 2.0.0 – 2026-10-09

Rewritten as an installable, model-agnostic package (`pip install -e .`,
command `snowtrack`) in response to the reviewers' request for a
reproducible "one-stop shop".

* General tracker: any number of soil layers (N ≥ 1), optional ice phase,
  any set of sources (`S_<NAME>` columns), three mixing conventions
  (`end`, `start`, `legacy`), unexplained water handled through the UNKNOWN
  source, per-step budget diagnostics; 30× faster than version 1.
* Adapters: VIC classic driver (`fluxes_*`), generic column mapping
  (YAML), netCDF (xarray), canonical tables; custom adapters via
  `register_adapter`.
* Pipeline: parallel runs, annual/monthly/overall summaries, run metadata.
* Routing of tracked components to a basin outlet with the Lohmann
  channel-routing scheme (wave velocities and diffusivity as calibration
  parameters), or with user-supplied distances / lags / weights.
* Metrics of the paper: fractions, T50 / centroid timing and trends,
  buffering index, late-season subsidy, apparent-lag diagnostics.
* Reproducibility: the archived per-cell outputs of Zhu et al. are
  reproduced to machine precision from raw VIC-glacier files
  (`--composition legacy`); continuous tests on Linux and Windows;
  Dockerfile; conda environment.
* Examples: minimal three-layer dual-phase cell (version 1 sample), four
  real VIC-glacier cells with archived outputs, a routed basin-outlet series
  of the study, a bring-your-own-model walk-through.

## 1.0.0 – 2026-06

Initial release: single-script three-layer dual-phase tracker with a
five-day sample (`snow_track/snow_track.py`). Reproduced in 2.0.0 with
`--composition start` (`examples/01_minimal_three_layer`).
