# Casebuilder

This directory owns source-side orchestration and prepared-case artifacts. It
does not decode GRIB or own runtime ingest logic inside the timestepper.

Current entry points:

- `check_external_toolchain.py`
- `prepare_case.py`
- `run_idealized_case.py`
- `toolchain_manifest.example.toml`

## Current Prepared-Case Contract

`prepare_case.py` now emits a concrete JSON artifact set for first-class source
adapters (`hrrr`, `rrfs`):

- `prepared_case_manifest.json`
- `analysis_state_stub.json`
- `boundary_cache_stub.json`
- `checkpoint_stub.json`
- `product_plan.json`

These artifacts mirror the strengthened ingest-side shape:

- analysis state:
  - atmosphere field group
  - surface field group
  - static-surface field group
- boundary cache:
  - monotonic snapshot offsets
  - atmosphere + surface groups
- checkpoint:
  - native restart/checkpoint stub contract only

The files are intentionally stubs:

- no real GRIB decode
- no source-data acquisition inside the runtime core
- no NetCDF writing

They exist to lock the casebuilder/runtime boundary now and give verification
tooling something concrete to validate.

## External Tooling Boundary

Preferred surrounding tooling remains:

- `rusbie` / `rustbie` for acquisition workflows
- `cfrust` for primary GRIB decode
- `ecrust` for low-level GRIB fallback

`prepare_case.py` reads the TOML toolchain manifest, records which external
repos are configured/present, and writes that availability into the prepared-case
manifest.

## Example

```powershell
python tools/casebuilder/prepare_case.py `
  --source hrrr `
  --domain-name central_plains_test `
  --forecast-hours 12 `
  --nx 240 --ny 180 --nz 70 `
  --dx 3000 --dy 3000 --z-top 22000
```

The output directory defaults to:

- `cases/prepared/<domain-name>_<source>/`

## Idealized Bridge

`run_idealized_case.py` remains the current concrete executable path into the
model. The new `product_plan.json` artifact explicitly records that
`idealized_summary_json` is the current verification/product bridge while real
source-driven map products are still being built out.
