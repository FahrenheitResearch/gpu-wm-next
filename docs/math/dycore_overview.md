# Dycore Overview

## Continuous equations

The MVP dycore is a compressible, nonhydrostatic finite-volume system on a
structured staggered grid.

## Discrete update

- slow tendencies are integrated with SSPRK3
- fast modes are delegated to a local split-explicit interface
- transport uses halo-aware face-flux differences
- thermal perturbations currently project onto vertical momentum through a
  buoyancy-only slow tendency
- horizontal pressure gradients now project onto `mom_u` and `mom_v` through
  face-centered dry pressure differences

## Invariants / admissibility

- constant-state preservation
- flux cancellation over interior faces
- closed-domain dry scalar conservation in periodic tests
- warm thermal perturbations must generate upward `w` response
- cold thermal perturbations must generate downward `w` response
- symmetric pressure perturbations must generate symmetric horizontal momentum
  response

## Assumptions for stability / consistency

- timestep respects the transport CFL assumptions for the chosen stencil
- halo exchange is correct before each stage requiring neighbor data

## Test mapping

- `tests/regression/test_constant_state.cpp`
- `tests/regression/test_scalar_mass_conservation.cpp`
- `tests/regression/test_hydrostatic_rest.cpp`
- `tests/regression/test_buoyancy_response.cpp`
- `tests/regression/test_horizontal_pressure_response.cpp`
- `tests/property/test_dry_momentum_virtual_rank_equivalence.cpp`
- `tests/unit/test_idealized_cases.cpp`
