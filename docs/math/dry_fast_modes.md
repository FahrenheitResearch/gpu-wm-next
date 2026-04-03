# Dry Fast Modes

## Continuous equations

The MVP fast path advances the dry compressible acoustic subsystem in
perturbation form:

- `d rho_d / dt = -div(m)`
- `d m_h / dt = -grad_h(p)`
- `d m_w / dt = -d p' / dz - rho' g`

where `m = (rho_d u, rho_d v, rho_d w)`, `p` is reconstructed from the dry
EOS using `rho_d` and `rho_d theta_m`, and the fast subcycle holds
`rho_d theta_m` fixed.

## Discrete update

- Use `dt_fast = dt / fast_substeps`.
- Reconstruct pressure each substep from the canonical conserved state.
- Update `mom_u` and `mom_v` from face-centered pressure differences.
- Update `mom_w` from a perturbation-form vertical pressure-gradient plus
  density-perturbation gravity term using:
  - horizontally averaged level reference profiles in the exact flat limit
  - column-local hydrostatic reference profiles when terrain-aware geometry is active
- Exchange face halos.
- Update `rho_d` from the divergence of the updated face momenta.

This is a local forward-backward acoustic subcycle. It does not yet include
nonlinear momentum transport, Coriolis, terrain metrics, or moist coupling.

## Invariants / admissibility conditions

- uniform state remains unchanged
- closed periodic domain conserves dry mass through the fast update
- hydrostatic-rest remains near rest under the perturbation-form vertical branch
- 1x1 and 2x2 virtual-rank layouts are equivalent within tolerance

## Assumptions for stability / consistency

- `fast_substeps >= 1`
- `dt_fast` is small enough for the explicit local acoustic update
- scalar and face halo exchange are correct before operators that need neighbor
  values
- the flat-limit reference reduces exactly to the prior 1-D profile behavior
- the terrain-aware reference is rebuilt from the local column state and metricized
  heights so terrain hydrostatic-rest remains near rest

## Test mapping

- `tests/unit/test_dry_fast_modes.cpp`
- `tests/property/test_dry_fastmode_virtual_rank_equivalence.cpp`
- `tests/regression/test_dry_constant_state_bundle.cpp`
- `tests/regression/test_hydrostatic_rest.cpp`
- `tests/regression/test_acoustic_pulse.cpp`
- `tests/regression/test_horizontal_pressure_response.cpp`
