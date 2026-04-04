# Warm-Rain Milestone

## Continuous equations

The first warm-rain milestone keeps moisture in the registry-managed tracer
state and adds local source terms for:

- water vapor: `rho_d q_v`
- cloud water: `rho_d q_c`
- rain water: `rho_d q_r`

The bounded source-term model is:

- `d/dt (rho_d q_v) = S_v`
- `d/dt (rho_d q_c) = S_c`
- `d/dt (rho_d q_r) = S_r`
- `S_v + S_c + S_r = 0`

Latent heating feeds back through `rho_d theta_m` during condensation and
evaporation. This checkpoint also adds bounded rain fallout to a separate
surface-accumulation reservoir, but it still avoids collection families beyond
the local warm-rain cut or any pressure-based moist state promotion.

## Discrete update

At the current checkpoint:

1. the dry core advances `rho_d`, `rho_d theta_m`, and face momentum
2. passive tracer transport advances the warm-rain tracer masses with the dry
   mass fluxes
3. a local warm-rain microphysics kernel updates `rho_d q_v`, `rho_d q_c`,
   `rho_d q_r`, and optionally `rho_d theta_m`
4. a bounded metric-aware fallout pass moves rain mass downward and deposits
   bottom-layer outflow into accumulated surface precipitation
5. the prepared-case runtime summary emits both dry diagnostics and tracer-side
   moisture totals/ranges, including accumulated precipitation closure

## Invariants / admissibility

- warm-rain tracers remain registry-managed and do not widen `DryState`
- total water mass (`q_v + q_c + q_r`) is conserved by the closed-box local
  warm-rain source update
- atmospheric water plus accumulated surface precipitation is conserved by the
  warm-rain-plus-fallout update
- tracer masses remain nonnegative
- runtime summary sidecars expose:
  - vapor, cloud, rain, condensed, and total water masses
  - accumulated/mean/max surface precipitation metrics
  - total surface cell count, precipitating cell count/fraction, and wet-cell
    mean accumulated precipitation
  - per-tracer total mass
  - per-tracer mixing-ratio min/max
- prepared/source-run summaries remain JSON-first artifacts, not restart truth

## Assumptions for stability / consistency

- the milestone is local and explicit
- fallout is first-order and uses the current hybrid-height cell thickness,
  not a placeholder representative layer depth
- no ice species, mixed phase, or aerosol activation
- no direct surface, PBL, or radiation coupling in the same patch
- source-run verification still treats the summary sidecar as a proof/report
  artifact rather than a restart contract
- bundle verification cross-checks the final summary fallout metrics against the
  emitted runtime plan-view accumulation/rate products when those artifacts are
  present

## Test mapping

- `tests/unit/test_warm_rain_microphysics.cpp`
- `tests/regression/test_warm_rain_closed_box.cpp`
- `tests/unit/test_runtime_summary.cpp`
- `tests/regression/test_warm_rain_summary_closure.cpp`
- `tools/verify/test_source_run_bundle.py`
- `tools/verify/run_verification.py`
