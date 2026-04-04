# Surface-Layer Closure

## Continuous equations

The shared surface-layer closure provides the neutral-limit baseline for
near-surface diagnostics and future surface-flux evaluation. In the neutral
limit, wind and scalar transfer follow logarithmic similarity relations with
separate momentum and heat roughness lengths:

- momentum interpolation uses `z0m`
- scalar interpolation uses `z0h`
- transfer coefficients are built from the corresponding log-profile factors

The closure consumes normalized physical forcing only:

- reference height and target heights
- reference wind, potential temperature, and humidity
- surface skin temperature and surface humidity
- surface pressure
- momentum and heat roughness lengths

## Discrete update

The current MVP computes one explicit neutral-limit diagnostic bundle:

- `u10`, `v10`
- `t2`, `q2`, `rh2`
- `cm`, `ch`, `cq`
- `ustar`
- reference and target wind-speed magnitudes

`screen_obsops` is a thin compatibility wrapper over this closure. The wrapper
must not own a second implementation of the near-surface math.

## Invariants / admissibility

- `z_ref`, `z_target_temp`, and `z_target_wind` must be positive
- `psfc` must remain physical
- effective `z0m` and `z0h` are positive after roughness clamping
- scalar diagnostics remain bounded between the surface and reference values in
  the neutral baseline
- with fixed temperature and pressure, increasing humidity inputs should not
  decrease neutral-limit `q2` or `rh2`
- `rh2` is clamped to `[0, 100]`
- wrapper and closure outputs must match exactly within floating-point
  tolerance

## Assumptions for stability / consistency

- this milestone implements the neutral-limit branch only
- future Monin-Obukhov stability corrections will extend this closure rather
  than replace it
- future surface-flux evaluation should reuse the same transfer-path outputs
  instead of duplicating similarity logic elsewhere

## Test mapping

- `tests/unit/test_surface_layer_closure.cpp`
- `tests/unit/test_surface_obsops.cpp`
- `tests/property/test_surface_tile_permutation.cpp`
- `tests/regression/test_surface_layer_neutral_reference.cpp`
