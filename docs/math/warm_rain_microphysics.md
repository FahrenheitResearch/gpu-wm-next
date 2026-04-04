# Warm-Rain Microphysics

## Scope

This is the first bounded warm-rain lane for `gpu-wm-next`. It is intentionally
local and explicit:

- prognostic water species remain tracer-managed masses on dry mass
- the operator updates:
  - `specific_humidity`
  - `cloud_water_mixing_ratio`
  - `rain_water_mixing_ratio`
- no ice
- no PBL/surface/radiation coupling
- simple rain fallout to a surface-accumulation sidecar

## Continuous model

For dry density `rho_d`, potential temperature `theta`, and liquid-water tracers
`q_v`, `q_c`, `q_r`:

- saturation adjustment:
  - if `q_v > q_sat(T, p)`, move excess vapor to cloud water
  - if `q_v < q_sat(T, p)` and `q_c > 0`, evaporate cloud water back to vapor
- autoconversion:
  - convert cloud water above a threshold to rain water
- optional rain evaporation:
  - if `q_v < q_sat(T, p)` and `q_r > 0`, evaporate a bounded amount of rain
    water back to vapor
- simple rain fallout:
  - move a bounded fraction of `q_r` downward each step
  - deposit bottom-layer outflow into accumulated surface precipitation

This first cut conserves atmospheric water locally when fallout is disabled:

`q_t = q_v + q_c + q_r`

With fallout enabled, the conserved quantity is atmospheric water plus the
surface-accumulation reservoir:

`q_t^* = q_v + q_c + q_r + P_s`

Latent heating is coupled only through vapor changes:

`Delta T = -(L_v / c_p) Delta q_v`

with potential-temperature update performed through the local Exner function.

## Discrete update

For each cell over timestep `Delta t`:

1. reconstruct `p` and `T` from `rho_d` and `rho_d theta`
2. compute `q_sat(T, p)`
3. apply saturation adjustment:
   - `C = max(q_v - q_sat, 0)`
   - `E_c = min(q_c, max(q_sat - q_v, 0))`
   - `q_v <- q_v - C + E_c`
   - `q_c <- q_c + C - E_c`
4. apply autoconversion:
   - `A = min(q_c, k_auto max(q_c - q_auto, 0) Delta t)`
   - `q_c <- q_c - A`
   - `q_r <- q_r + A`
5. apply bounded rain evaporation:
   - `E_r = min(q_r, k_evap min(q_sat - q_v, q_r) Delta t)`
   - `q_r <- q_r - E_r`
   - `q_v <- q_v + E_r`
6. apply bounded fallout when terminal velocity is positive:
   - if terrain-aware metrics are supplied, compute per-cell
     `Delta z(i, j, k)` from the hybrid-height geometry
   - otherwise use the configured representative
     `sedimentation_layer_depth_m`
   - `f(i, j, k) = min(v_t Delta t / Delta z(i, j, k), 1)`
   - move rain-column mass per area downward one layer at a time
   - accumulate bottom-layer outflow in `P_s` with units `kg m^-2` / `mm`
7. if latent-heating coupling is enabled:
   - `Delta q_v = q_v^{n+1} - q_v^n`
   - `theta <- theta - (L_v / c_p) Delta q_v / Pi`

Tracer storage remains conservative in `rho_d q_x`.

## Invariants / admissibility

- nonnegativity of `q_v`, `q_c`, `q_r`
- local total-water conservation for the source operator
- column water plus surface accumulation conservation when fallout is enabled
- latent heating only from vapor-phase change
- no changes to dry density or momentum in this operator

## Assumptions / limitations

- uses dry pressure/temperature reconstruction
- fallout is first-order and bounded, but it now uses the actual hybrid-height
  cell thickness rather than a representative layer depth
- no accretion in this cut
- no supersaturation iteration; one explicit local adjustment only
- intended as the first moist-product lane, not a full operational warm-rain
  package

## Test mapping

- `tests/unit/test_warm_rain_microphysics.cpp`
  - supersaturated cell condenses vapor to cloud water
  - autoconversion shifts cloud water to rain water
  - local total-water conservation
  - latent heating raises `theta`
  - bounded representative-layer fallout preserves nonnegativity and column
    closure in the local overload
  - single-column rain fallout transfers rain water into the surface bucket
- `tests/regression/test_warm_rain_closed_box.cpp`
  - closed-box water conservation over multiple cells
  - positive cloud/rain production from a moist warm box
  - increased `rho_d theta` from condensation heating
  - atmospheric rain plus surface accumulation is conserved under fallout
- `tests/regression/test_warm_rain_summary_closure.cpp`
  - runtime summary sidecar preserves total water plus accumulated
    precipitation closure and nonnegative tracer ranges
