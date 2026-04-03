# Surface Fluxes

## Continuous equations

The intended architecture solves a surface energy and moisture balance using
surface state, soil/canopy state, and atmospheric forcing.

## Discrete update

The Stage 0 implementation only commits the data model and neutral-limit screen
operator contract. Full flux evolution is deferred to Stage 2.

## Invariants / admissibility

- soil moisture remains bounded
- canopy water remains bounded
- surface diagnostic operators remain explicit and separate from proxy fields

## Assumptions for stability / consistency

- a surface energy-balance solver will later share closures with the screen
  operators
- tile state is local to a horizontal cell and does not require horizontal halo
  exchange

## Test mapping

- current: `tests/unit/test_surface_obsops.cpp`
- future: energy residual, moisture boundedness, and flux monotonicity tests
