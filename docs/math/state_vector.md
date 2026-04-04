# State Vector

## Continuous equations

The MVP atmospheric state is evolved in conservative form for dry air, momentum,
moist entropy surrogate, and tracers:

- `∂t rho_d + div(rho_d u) = 0`
- `∂t (rho_d u) + div(F_m) = S_m`
- `∂t (rho_d theta_m) + div(F_theta) = S_theta`
- `∂t (rho_d q_k) + div(F_k) = S_k`

## Discrete update

Finite-volume cell averages are stored on a halo-padded structured grid. Face
flux differences update owned cells. Tracer masses are updated by the same
transport machinery as other conserved quantities.

For the first bounded moisture milestone, the canonical transported moisture
name is `specific_humidity`. The ingest/runtime bridge may still accept
`water_vapor_mixing_ratio` as a transitional source alias, but that alias is
not part of the runtime state contract and should not appear as the canonical
product name in tracer-aware outputs.

At the current checkpoint, the source-driven product bridge can already append
moist diagnostics from populated companion analysis data, and the first
registry-managed warm-rain tracer set is active in the runtime:

- `specific_humidity`
- `cloud_water_mixing_ratio`
- `rain_water_mixing_ratio`

That does not change the canonical runtime naming contract above.

## Invariants / admissibility

- dry-mass conservation on a closed domain
- constant-state preservation
- tracer non-negativity when a positivity limiter is enabled
- canonical tracer names remain stable across the runtime/product boundary
- warm-rain source terms conserve total water in the closed-box milestone cut

## Assumptions for stability / consistency

- consistent face metrics
- halo values reflect the active boundary contract
- tracer fluxes use the same geometry/sign conventions as dry mass
- first moist diagnostic products derived from tracer state are products, not
  canonical prognostic state
- warm-rain microphysics stays local and does not widen `DryState`

## Test mapping

- `tests/regression/test_constant_state.cpp`
- `tests/regression/test_scalar_mass_conservation.cpp`
- `tests/unit/test_tracer_registry.cpp`
- `tests/unit/test_warm_rain_microphysics.cpp`
- `tests/regression/test_warm_rain_closed_box.cpp`
- `tests/regression/test_warm_rain_summary_closure.cpp`
