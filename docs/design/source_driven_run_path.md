# Source-Driven Run Path

This note describes the prepared-case bridge from source data into runtime
artifacts and then into map products. It is intentionally a contract note, not
a runtime implementation spec.

## Goal

Make a prepared-case directory usable as the handoff point between source
decode, runtime execution, verification, and map rendering.

The bridge is:

1. `prepare_case.py` writes the prepared-case manifest and artifact contracts.
2. `populate_prepared_case.py` fills the analysis and boundary JSON artifacts
   from external source tooling.
3. The runtime driver consumes `analysis_state.json` and
   `boundary_cache.json`.
4. The runtime emits summary, plan-view, and map-manifest JSON artifacts.
5. When a populated companion `analysis_state.json` is present, the current
   plan-view writer may enrich the source-driven bundle with moisture
   diagnostics such as `specific_humidity`, `relative_humidity`, and
   `dewpoint`.
6. `tools/verify/` validates the bundle.

## Bundle Contents

The populated bundle is expected to carry these JSON artifacts:

- `prepared_case_manifest.json`
- `analysis_state.json`
- `boundary_cache.json`
- `summary.json`
- `plan_view.json`
- `map_manifest.json`

The verifier also still supports the stub-only prepared-case manifest path so
the contract can be checked before population.

## Verification Boundary

`tools/verify/run_verification.py` validates:

- prepared-case manifests
- populated analysis-state artifacts
- populated boundary-cache artifacts
- prepared-case/source-run bundle directories
- plan-view and map-manifest outputs

This keeps source-driven orchestration outside the timestepper while still
making the handoff concrete enough to test.

## Current Constraints

- source decode remains outside the runtime core
- runtime code only consumes populated canonical JSON artifacts
- verification is JSON-first and does not decode GRIB itself
- the source-run path should keep the dry-core map bridge available even before
  full moist physics lands

## Test Mapping

- `tools/verify/test_source_run_bundle.py`
- `tools/verify/run_verification.py`
- `tools/verify/render_plan_view_maps.py`
- `docs/math/prepared_case_runtime.md`
