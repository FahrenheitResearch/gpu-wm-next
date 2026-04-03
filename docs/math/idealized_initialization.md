# Idealized Initialization

## Continuous equations

Idealized dry cases prescribe analytic or semi-analytic background states plus
localized perturbations that map into the canonical dry conserved state.

## Discrete update

- build a rectilinear domain and vertical metric set
- initialize conserved dry mass and dry potential-temperature density
- initialize face momentum from prescribed background flow
- apply localized thermal perturbations through `rho_d * theta_m`
- keep perturbation support explicit and bounded

## Invariants / admissibility

- initialized dry density remains positive everywhere
- initialized `rho_d * theta_m` remains positive everywhere
- perturbation support is bounded to the configured ellipsoid
- flat background states reduce to the constant/hydrostatic base cases

## Assumptions for stability / consistency

- background density/theta values are physically positive
- perturbation amplitudes do not drive `theta_m` negative
- domain metrics are consistent with the declared horizontal/vertical extents

## Test mapping

- `tests/unit/test_idealized_domain_builder.cpp`
- `tests/unit/test_idealized_topography.cpp`
- `tests/unit/test_idealized_cases.cpp`
