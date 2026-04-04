# Casebuilder

This directory owns source-side orchestration and prepared-case artifacts. It
does not decode GRIB or own runtime ingest logic inside the timestepper.

Current entry points:

- `check_external_toolchain.py`
- `prepare_case.py`
- `populate_prepared_case.py`
- `run_prepared_case.py`
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

## Runtime Bridges

`run_idealized_case.py` remains the fastest local dry-core path into the model.
The prepared/source-driven path is now also executable:

1. `prepare_case.py` writes the prepared-case manifest and artifact contract.
2. `populate_prepared_case.py` fills the analysis and boundary JSON artifacts
   from the configured external toolchain.
3. `run_prepared_case.py` drives `gwm_prepared_case_driver`, renders plan-view
   maps, and runs verification on the resulting bundle.

Both wrappers can request a plan-view field slice from the driver and render
actual map images in one step.

Example:

```powershell
python tools/casebuilder/run_idealized_case.py `
  cases/idealized/warm_bubble.yaml `
  --steps 2 `
  --dt 2.0 `
  --fast-substeps 2 `
  --plan-view-level 6 `
  --render-maps
```

That produces:

- `<case>.summary.json`
- `<case>.plan_view.json`
- `<case>.plan_view_maps/map_manifest.json`
- one rendered image per emitted field

Prepared/source-driven example:

```powershell
python tools/casebuilder/prepare_case.py `
  --source hrrr `
  --domain-name central_plains_smoke `
  --cycle-time-utc 2026-04-04T00:00:00Z `
  --forecast-hours 1 `
  --nx 24 --ny 24 `
  --center-lat 39.0 --center-lon -97.0 `
  --pressure-levels-hpa 1000,850,700

python tools/casebuilder/run_prepared_case.py `
  --prepared-case cases/prepared/central_plains_smoke_hrrr/prepared_case_manifest.json `
  --populate `
  --steps 1 `
  --dt 2.0 `
  --fast-substeps 2 `
  --plan-view-level 1 `
  --render-maps
```

That produces a verified source-run bundle containing:

- `prepared_case_manifest.json`
- `analysis_state.json`
- `boundary_cache.json`
- `summary.json`
- `plan_view.json`
- `map_manifest.json`
