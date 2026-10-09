# FAQ

**My model has no explicit infiltration variable.** Leave `INFIL` out; it is
taken as `max(sum(S_*) - RUNOFF_SURF, 0)`.

**My model has one soil moisture value, no layers.** That is the single-store
configuration (`W1_LIQ`, `Q1_BASE`, optionally `E1_TR`).

**Which `--composition` should I use?** `end` (default) for new work;
`legacy` to reproduce the archived outputs of Zhu et al.; `start` to
reproduce version 1 of this repository. They differ only in the timing of
mixing within a step (`docs/method.md`, section 3).

**How long must the spin-up be?** Watch `FW_UNKNOWN` and
`RUNOFF_UNKNOWN_TOTAL`: the UNKNOWN share decays as the stores turn over
(years to decades for deep stores). Either track from a date well before the
analysis period and drop the output with `--spinup-years N`, or drop the
input years before tracking with `--skip-input-years N` (archived convention).

**The closure error is 1e-8 instead of 0.** Stores smaller than `epsilon`
(1e-8 mm) are treated as empty; this bounds the closure error by `epsilon`.

**Runoff of some days is attributed to UNKNOWN even after spin-up.** Those
are days on which direct runoff exceeded the water that reached the ground
(or infiltration exceeded it); the excess cannot be attributed to a source.
`--strict-surface-closure` turns this into an error if you prefer.

**Can I track sources other than rain / snow / glacier?** Yes: any column
`S_<NAME>` is a source (e.g. `S_PERMAFROST`, `S_IRRIGATION`).

**Sub-daily time steps?** Provide a `DATE` column (or `HOUR`); all fluxes are
per step. Routing kernels are then in steps of that length.

**Which routing velocities should I use?** They are calibration parameters:
calibrate the surface and subsurface wave velocities (and, if needed, the
diffusivity) of `snowtrack route` against gauged discharge together with your
model, or supply the lags of your own routing model through `--cells`.

**Large domains?** `--workers N` parallelises over cells; `--columns compact
--float-format %.8g --gzip` keeps the output small (about 3 s per 60-year
daily cell and core).

**Is this a water-age or tracer model?** No. Compositions are per store and
step (complete mixing); the lag diagnostics are effective storage-release
time scales of the host model, not isotope-equivalent water ages.
