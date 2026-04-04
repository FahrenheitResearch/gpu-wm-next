# Moisture Tracer And Products

## Continuous equations

The first bounded moisture milestone promotes water vapor to the first
registry-managed tracer while keeping the dry-core stepping architecture
otherwise unchanged. In conservative form, the intended transport equation is:

- `d/dt (rho_d q_v) + div(F_qv) = S_qv`

For the passive first cut:

- `S_qv = 0`

At the current checkpoint, source-driven moist plan-view products can already
be derived diagnostically from populated companion analysis data together with
the local thermodynamic state. The longer-term intended runtime/product fields
are:

- canonical moisture field: `specific_humidity`
- derived diagnostics: `relative_humidity`, `dewpoint`

## Discrete update

The runtime/product contract for this milestone is:

1. source-driven artifacts already provide populated humidity inputs under the
   canonical `specific_humidity` name, with a transitional ingest alias still
   accepted at the bridge edge
2. the tracer registry carries the canonical runtime tracer name
3. the runtime transport path should treat the first moisture tracer like any
   other conserved tracer mass once tracer stepping is threaded into the core
4. plan-view/map products extend the existing `gwm-next-plan-view/v1` bundle
   with moisture diagnostics at the same `nx * ny` slice shape and row-major
   storage semantics as current dry fields

## Invariants / admissibility

- `specific_humidity` is the canonical runtime/product moisture name
- `water_vapor_mixing_ratio` remains a source-ingest alias only
- moisture-tracer storage preserves declared units and positivity metadata
- moist plan-view fields must satisfy `len(values) == nx * ny`
- moist diagnostic outputs must remain finite for admissible thermodynamic
  inputs
- the moisture extension must not create a second canonical state container for
  diagnostics

## Assumptions for stability / consistency

- the first moisture milestone is passive first:
  - no latent-heating coupling
  - no condensation/evaporation source terms
  - no microphysics in the same patch
- current source-driven moist diagnostics come from populated companion
  analysis-state inputs; full runtime tracer transport is still follow-on work
- derived `relative_humidity` and `dewpoint` are diagnostic products, not
  restart truth
- near-surface humidity diagnostics continue to come through the shared
  surface-layer closure / screen-obsops path

## Test mapping

- `tests/unit/test_tracer_registry.cpp`
- `tests/unit/test_surface_obsops.cpp`
- `tests/unit/test_ingest_contracts.cpp`
- `tests/unit/test_runtime_case.cpp`
- future tracer-transport and moist plan-view tests should extend the existing
  runtime/product suites once the runtime stepping grows beyond the current
  source-driven enrichment path
