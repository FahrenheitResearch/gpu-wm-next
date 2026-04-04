# Moisture Tracer And Products

## Continuous equations

The first bounded moisture milestone promotes water vapor to the first
registry-managed tracer and then extends that runtime contract to a warm-rain
tracer set while keeping the dry-core ownership architecture unchanged. In
conservative form, the transport/source equations are:

- `d/dt (rho_d q_v) + div(F_qv) = S_qv`
- `d/dt (rho_d q_c) + div(F_qc) = S_qc`
- `d/dt (rho_d q_r) + div(F_qr) = S_qr`

For the passive first cut:

- `S_qv = 0`

For the current warm-rain milestone:

- `S_qv + S_qc + S_qr = 0`
- latent heating may update `rho_d theta_m`, but no new prognostic pressure or
  moist-state container is introduced

At the current checkpoint, source-driven moist plan-view products can already
be derived diagnostically from populated companion analysis data together with
the local thermodynamic state. The longer-term intended runtime/product fields
are:

- canonical moisture field: `specific_humidity`
- warm-rain tracer products:
  - `cloud_water_mixing_ratio`
  - `rain_water_mixing_ratio`
  - `total_condensate`
  - `column_rain_water`
- derived diagnostics:
  - `relative_humidity`
  - `dewpoint`
  - `synthetic_reflectivity`

## Discrete update

The runtime/product contract for this milestone is:

1. source-driven artifacts already provide populated humidity inputs under the
   canonical `specific_humidity` name, with a transitional ingest alias still
   accepted at the bridge edge
2. the tracer registry carries the canonical runtime tracer name
3. the runtime transport path treats the warm-rain tracer set like any other
   conserved tracer masses once tracer stepping is threaded into the core
4. plan-view/map products extend the existing `gwm-next-plan-view/v1` bundle
   with moisture diagnostics at the same `nx * ny` slice shape and row-major
   storage semantics as current dry fields
5. warm-rain plan-view products are derived from runtime tracer masses, not
   clobbered by prepared-case enrichment sidecars
6. the prepared-case summary sidecar reports tracer totals plus total/condensed
   water masses for milestone-proof closure

## Invariants / admissibility

- `specific_humidity` is the canonical runtime/product moisture name
- `water_vapor_mixing_ratio` remains a source-ingest alias only
- moisture-tracer storage preserves declared units and positivity metadata
- the warm-rain tracer set is:
  - `specific_humidity`
  - `cloud_water_mixing_ratio`
  - `rain_water_mixing_ratio`
- moist plan-view fields must satisfy `len(values) == nx * ny`
- moist diagnostic outputs must remain finite for admissible thermodynamic
  inputs
- `synthetic_reflectivity` is emitted only when rain water is present in the
  runtime tracer contract; there is no fake source-side substitute
- the moisture extension must not create a second canonical state container for
  diagnostics

## Assumptions for stability / consistency

- current source-driven moist diagnostics still include derived plan-view
  products from populated companion analysis-state inputs
- warm-rain remains bounded:
  - no sedimentation/fallout yet
  - no mixed-phase ice species
  - no surface/PBL/radiation coupling in the same patch
- derived `relative_humidity` and `dewpoint` are diagnostic products, not
  restart truth
- near-surface humidity diagnostics continue to come through the shared
  surface-layer closure / screen-obsops path

## Test mapping

- `tests/unit/test_tracer_registry.cpp`
- `tests/unit/test_warm_rain_microphysics.cpp`
- `tests/unit/test_runtime_summary.cpp`
- `tests/regression/test_warm_rain_closed_box.cpp`
- `tests/regression/test_warm_rain_summary_closure.cpp`
- `tests/unit/test_surface_obsops.cpp`
- `tests/unit/test_ingest_contracts.cpp`
- `tests/unit/test_runtime_case.cpp`
- `tools/verify/test_source_run_bundle.py`
