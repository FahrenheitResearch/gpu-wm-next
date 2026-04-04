# Verification

This directory owns successor-model validation and product-side JSON plumbing.
It validates contracts and current idealized outputs; it does not replace the
future diagnostics/plotting stack.

Current entry points:

- `run_verification.py`
- `render_plan_view_maps.py`

## Supported Inputs

`run_verification.py` now understands:

- `prepared_case_manifest.json`
- populated `analysis_state.json`
- populated `boundary_cache.json`
- `boundary_cache_stub.json`
- `checkpoint_stub.json`
- `product_plan.json`
- plan-view JSON from `gwm_idealized_driver`
- `map_manifest.json`
- idealized driver summary JSON from `run_idealized_case.py`
- source-run bundle directories containing:
  - `prepared_case_manifest.json`
  - `analysis_state.json`
  - `boundary_cache.json`
  - optional `summary.json`
  - optional `plan_view.json`
  - optional `map_manifest.json`

Each verification run emits a concrete report JSON next to the input by default.

## Current Verification Scope

Prepared-case side:

- schema/version checks
- required field-group presence
- monotonic boundary offsets
- referenced artifact existence
- external toolchain availability summary
- populated analysis/boundary artifact checks
- source-run bundle checks for prepared-case directories

Idealized side:

- summary-shape checks
- finite mass/theta totals
- finite `w_face` extrema
- derived drift/span metrics
- explicit bridge to current product support
- plan-view bundle shape checks
- rendered image-path checks from `map_manifest.json`

## Examples

Verify a prepared-case manifest:

```powershell
python tools/verify/run_verification.py `
  --input cases/prepared/central_plains_test_hrrr/prepared_case_manifest.json
```

Verify an idealized dry-core summary:

```powershell
python tools/verify/run_verification.py `
  --input cases/idealized/warm_bubble.summary.json
```

Verify a rendered map manifest:

```powershell
python tools/verify/run_verification.py `
  --input cases/idealized/warm_bubble.plan_view_maps/map_manifest.json
```

Verify a prepared-case/source-run bundle directory:

```powershell
python tools/verify/run_verification.py `
  --input cases/prepared/<populated-case-directory> `
  --kind source_run_bundle
```

Render plan-view maps directly:

```powershell
python tools/verify/render_plan_view_maps.py `
  --input cases/idealized/warm_bubble.plan_view.json
```

Smoke-test the source-run verification path:

```powershell
python tools/verify/test_source_run_bundle.py
```

## External Tooling Boundary

These scripts are intentionally JSON-first and orchestration-oriented:

- no GRIB decode
- no runtime-core field mutation
- plotting here is limited to the current plan-view JSON product

They are designed to hand off later to:

- `wrf-rust` for diagnostics / verification
- `wrf-rust-plots` for plots / panels / map products
- Rust thermo/severe tools for derived diagnostics
