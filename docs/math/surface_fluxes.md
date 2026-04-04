# Surface Fluxes

## Continuous equations

The intended architecture solves a surface energy and moisture balance using
surface state, soil/canopy state, and atmospheric forcing.

## Discrete update

The current increment commits:

- mutable `SurfaceState`
- immutable `SurfaceStaticProperties`
- a neutral-limit shared surface-layer closure
- tile aggregation through `surface_layer_exchange`

Full flux evolution and atmosphere coupling are still deferred.

## Invariants / admissibility

- soil moisture remains bounded
- canopy water remains bounded
- surface diagnostic operators remain explicit and separate from proxy fields

## Assumptions for stability / consistency

- a surface energy-balance solver will later share closures with the screen
  operators
- tile state is local to a horizontal cell and does not require horizontal halo
  exchange
- the current closure milestone is limited to neutral diagnostics and
  transfer-path outputs; it does not yet couple flux tendencies into the
  atmospheric timestep

## Test mapping

- `tests/unit/test_surface_state.cpp`
- `tests/unit/test_surface_static_properties.cpp`
- `tests/unit/test_surface_layer_closure.cpp`
- `tests/unit/test_surface_obsops.cpp`
- `tests/property/test_surface_tile_permutation.cpp`
- `tests/regression/test_surface_layer_neutral_reference.cpp`
- future: energy residual, moisture boundedness, and flux monotonicity tests
