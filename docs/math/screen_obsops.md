# Screen Observation Operators

## Continuous equations

The screen operators are based on Monin-Obukhov similarity, mapping lowest-model
level and surface state into 2 m / 10 m quantities.

## Discrete update

The current screen operators are a thin wrapper over the shared
`surface_layer_closure` neutral-limit path. They do not own separate similarity
math.

## Invariants / admissibility

- screen wind should decrease monotonically as target height decreases, all else
  equal
- scalar interpolation should remain bounded between surface and lowest-model
  level values in the neutral-limit baseline
- with fixed temperature and pressure, increasing the reference/surface humidity
  inputs should not decrease neutral-limit `q2` or `rh2`
- wrapper outputs must match the shared closure outputs exactly within floating
  point tolerance

## Assumptions for stability / consistency

- target height is above effective roughness length
- surface and lowest-model-level states are physically ordered
- stability corrections remain deferred; the neutral-limit branch is the current
  analytic baseline

## Test mapping

- `tests/unit/test_surface_obsops.cpp`
- `tests/unit/test_surface_layer_closure.cpp`
- `tests/property/test_surface_tile_permutation.cpp`
- `tests/regression/test_surface_layer_neutral_reference.cpp`
- future stable/unstable Monin-Obukhov limit tests
