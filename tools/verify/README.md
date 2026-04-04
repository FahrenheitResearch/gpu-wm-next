# Verification

This directory owns successor-model validation and product-side JSON plumbing.
It validates contracts and current idealized outputs; it does not replace the
future diagnostics/plotting stack.

Current entry point:

- `run_verification.py`

## Supported Inputs

`run_verification.py` now understands:

- `prepared_case_manifest.json`
- `boundary_cache_stub.json`
- `checkpoint_stub.json`
- `product_plan.json`
- idealized driver summary JSON from `run_idealized_case.py`

Each verification run emits a concrete report JSON next to the input by default.

## Current Verification Scope

Prepared-case side:

- schema/version checks
- required field-group presence
- monotonic boundary offsets
- referenced artifact existence
- external toolchain availability summary

Idealized side:

- summary-shape checks
- finite mass/theta totals
- finite `w_face` extrema
- derived drift/span metrics
- explicit bridge to current product support

## Example

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

## External Tooling Boundary

This script is intentionally JSON-first and orchestration-oriented:

- no GRIB decode
- no runtime-core field mutation
- no plotting implementation here

It is designed to hand off later to:

- `wrf-rust` for diagnostics / verification
- `wrf-rust-plots` for plots / panels / map products
- Rust thermo/severe tools for derived diagnostics
