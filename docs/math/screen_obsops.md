# Screen Observation Operators

## Continuous equations

The screen operators are based on Monin-Obukhov similarity, mapping lowest-model
level and surface state into 2 m / 10 m quantities.

## Discrete update

The MVP scaffold uses a neutral-limit operator with explicit roughness-length
dependence as the first analytic baseline.

## Invariants / admissibility

- screen wind should decrease monotonically as target height decreases, all else
  equal
- scalar interpolation should remain bounded between surface and lowest-model
  level values in the neutral-limit baseline

## Assumptions for stability / consistency

- target height is above effective roughness length
- surface and lowest-model-level states are physically ordered

## Test mapping

- `tests/unit/test_surface_obsops.cpp`
- future stable/unstable Monin-Obukhov limit tests
