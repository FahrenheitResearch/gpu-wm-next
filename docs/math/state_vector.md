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

## Invariants / admissibility

- dry-mass conservation on a closed domain
- constant-state preservation
- tracer non-negativity when a positivity limiter is enabled

## Assumptions for stability / consistency

- consistent face metrics
- halo values reflect the active boundary contract
- tracer fluxes use the same geometry/sign conventions as dry mass

## Test mapping

- `tests/regression/test_constant_state.cpp`
- `tests/regression/test_scalar_mass_conservation.cpp`
- future hydrostatic-rest and tracer-positivity tests
