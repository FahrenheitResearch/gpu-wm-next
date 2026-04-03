# Dry Pressure Gradient

## Continuous equations

For the dry momentum MVP, horizontal momentum tendencies include a pressure
gradient term of the form:

- `d(rho_d u)/dt = -dp/dx`
- `d(rho_d v)/dt = -dp/dy`

with dry pressure diagnosed from the canonical dry state through the dry ideal
gas relation written in terms of `rho_d` and `theta_m`.

## Discrete update

- diagnose cell-centered dry pressure from `rho_d` and `rho_d * theta_m`
- use halo-exchanged pressure on the local cell-centered stencil
- evaluate face-centered pressure differences on the C-grid
- add x-face tendencies to `mom_u`
- add y-face tendencies to `mom_v`
- keep open outer-domain faces zeroed for this MVP when no lateral neighbor
  exists

## Invariants / admissibility

- uniform dry pressure produces zero horizontal pressure-gradient tendency
- symmetric centered thermal perturbations produce antisymmetric horizontal
  momentum tendencies
- hydrostatic rest with horizontally uniform state keeps `mom_u` and `mom_v`
  near zero
- serial and virtual-rank decompositions produce equivalent pressure-gradient
  updates

## Assumptions for stability / consistency

- the dry EOS helper is applied consistently from the canonical state only
- cell-centered pressure halos are valid before face tendencies are evaluated
- x/y face halos follow the decomposition contract and are exchanged without
  host-side stitching
- this MVP operator excludes Coriolis, map-factor terms, nonlinear momentum
  advection, and split-explicit acoustic handling

## Test mapping

- `tests/unit/test_face_halo_exchange.cpp`
- `tests/unit/test_dry_pressure_gradient.cpp`
- `tests/property/test_dry_momentum_virtual_rank_equivalence.cpp`
- `tests/regression/test_hydrostatic_rest.cpp`
- `tests/regression/test_horizontal_pressure_response.cpp`
