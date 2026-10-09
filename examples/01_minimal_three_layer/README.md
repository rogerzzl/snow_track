# Example 01 – minimal three-layer, dual-phase cell

A five-day hand-made cell in the canonical schema (`input/grid_001.csv`): three
soil layers, liquid **and** ice storages, per-layer percolation / baseflow /
transpiration and freeze–thaw fluxes, three sources (rain, snow, glacier).
It is the smallest input that exercises every term of the scheme in
`docs/method.md`.

```bash
snowtrack run -i examples/01_minimal_three_layer/input -o examples/01_minimal_three_layer/output --columns full
```

`expected/tracked_grid_001_v1_start.csv` is the output of version 1 of this
repository, which used the explicit start-of-step composition. It is reproduced
to machine precision with

```bash
snowtrack run -i examples/01_minimal_three_layer/input -o out_v1 --columns full --composition start --strict-surface-closure
```

(`tests/test_tracker.py::test_v1_three_layer_sample_is_reproduced_with_start_composition`).
The default composition (`end`, mix-then-withdraw) gives slightly different
numbers; see `docs/method.md` for the three conventions.
