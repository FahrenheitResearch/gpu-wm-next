# Dry Momentum Flux

## Continuous equations

The slow dry-momentum transport extension advances the conservative momentum
components in flux form on the existing Cartesian C-grid:

- `d m_u / dt = -div(F_u)`
- `d m_v / dt = -div(F_v)`
- `d m_w / dt = -div(F_w)`

where each `m_*` is the corresponding face-centered momentum component, and the
transport fluxes are formed from the current momentum state and cell-centered
advecting velocities reconstructed from the canonical conserved variables.

## Discrete update

- Reconstruct haloed cell-centered velocities from `rho_d` and face-centered
  momentum fields.
- For each face-oriented momentum field, compute a first-order upwind
  flux-divergence tendency on its native control volume.
- Add that tendency into the existing slow operator alongside:
  - continuity/theta transport
  - horizontal pressure-gradient source
  - vertical buoyancy source

This first pass is intentionally monotone and low-order. It does not yet
include terrain metrics, Coriolis, damping, diffusion, or higher-order
reconstruction.

## Invariants / admissibility conditions

- constant state remains constant
- periodic dry-momentum transport remains near-conservative
- 1x1 and 2x2 virtual-rank layouts remain equivalent within tolerance
- hydrostatic rest remains near rest once transport is added
- short density-current evolution remains finite and symmetry-preserving

## Assumptions for stability / consistency

- timestep satisfies the current explicit transport CFL assumptions
- scalar and face halos are exchanged before constructing the momentum-flux
  tendencies
- cell-centered velocity scratch is derived only from the canonical conserved
  state and remains ephemeral scratch, not new model truth

## Test mapping

- `tests/unit/test_dry_momentum_flux.cpp`
- `tests/property/test_dry_momentum_virtual_rank_equivalence.cpp`
- `tests/regression/test_constant_state.cpp`
- `tests/regression/test_hydrostatic_rest.cpp`
- `tests/regression/test_buoyancy_response.cpp`
- `tests/regression/test_horizontal_pressure_response.cpp`
- `tests/regression/test_density_current_evolution.cpp`
