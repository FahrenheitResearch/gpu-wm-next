# Plan-View Products

## Continuous quantities

The first rendered map path is intentionally diagnostic, not prognostic. It
extracts plan-view fields from the current dry-core state for visualization and
verification:

- `rho_d`
- `theta_m = rho_theta_m / rho_d`
- cell-centered velocity estimates from C-grid face momenta
- derived horizontal wind speed
- terrain height
- physical slice height

## Discrete extraction

The current implementation:

1. gathers distributed dry-state fields into a global plan-view bundle using the
   existing virtual-rank layout contract
2. derives cell-centered velocity estimates from adjacent face momenta divided
   by local dry density
3. writes a JSON bundle with explicit `nx`, `ny`, row-major storage, units, and
   a list of named fields
4. feeds that bundle into Python-side map rendering and manifest generation

This is deliberately not a restart/checkpoint format and deliberately not a
general I/O backend. It is the first model-to-map bridge only.

## Prepared-case moist enrichment

The current source-driven product bridge can enrich the dry plan-view bundle
from a populated companion `analysis_state.json` while keeping the same bundle
semantics. The current moisture-aware diagnostics are:

- `specific_humidity` as the canonical transported moisture field
- `relative_humidity` derived from local thermodynamic state
- `dewpoint` derived from local thermodynamic state

The bridge should not expose source aliases such as
`water_vapor_mixing_ratio` as plan-view product names. Moist products remain
diagnostic outputs layered on top of the runtime state; they are not a second
canonical state vector and they are not a restart schema.

## Invariants / admissibility

- output schema version must be explicit
- each field must satisfy `len(values) == nx * ny`
- slice index must remain within `[0, nz - 1]`
- gathered plan-view fields must preserve constant-state values
- rendered map manifests must reference existing image artifacts
- moist extensions must preserve the same row-major field sizing and slice
  semantics as the existing dry bundle

## Assumptions for stability / consistency

- current plan-view products are for the final dry state of a run
- moisture enrichment extends the same plan-view schema rather than inventing a
  second map-product container
- current cell-centered velocity diagnostics are visualization-friendly derived
  products, not canonical prognostic state
- current rendering is idealized-case focused, but the same plan-view schema is
  now used as the source-driven prepared-case bridge into map products
- real-data map parity still depends on the source-run runtime path, not just
  on the renderer

## Test mapping

- `tests/unit/test_plan_view_output.cpp`
- `tests/unit/test_tracer_registry.cpp`
- `tools/verify/run_verification.py`:
  - `plan_view_bundle`
  - `map_manifest`
- source-driven bundle verification path:
  - `tools/verify/test_source_run_bundle.py`
  - `tools/verify/run_verification.py` on prepared-case/source-run directories
- end-to-end smoke path through:
  - `tools/casebuilder/run_idealized_case.py`
  - `tools/verify/render_plan_view_maps.py`
- future moisture/tracer product checks should extend the existing plan-view
  tests rather than fork a separate map schema
