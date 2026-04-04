# External Toolchain Integration

This repository intentionally separates:

- runtime core: C++20 + CUDA + MPI
- surrounding ecosystem: Rust/Python tools for acquisition, ingest, diagnostics,
  plotting, and verification

## High-Level Policy

- The core model owns:
  - state
  - dycore
  - physics
  - surface/land model
  - decomposition / halos
  - checkpoint / restart state
- External tooling owns:
  - pulling source data
  - high-level GRIB decoding
  - low-level GRIB mutation / sidecar workflows
  - postprocessing / severe-weather diagnostics
  - product rendering

## Preferred External Repos

### Case acquisition / ingest

- `rusbie` / `rustbie`: fast source-data acquisition
- `cfrust`: primary high-level GRIB ingest path
- `ecrust`: low-level GRIB/ecCodes-style fallback

### Verification / products

- `wrf-rust`: diagnostics / verification
- `wrf-rust-plots`: plotting / product generation
- Rust MetPy replacement repo: generic thermodynamic / severe-weather math

## Intended Repo Boundary

- `tools/casebuilder/` prepares target-area cases and boundary caches by
  orchestrating external data tooling.
- `tools/verify/` runs successor-model validation and user-facing product
  generation by orchestrating diagnostics/plotting tools.
- `src/ingest/` defines canonical interfaces only. It does not absorb those
  external repos into the timestepper.

## Current Artifact Hand-Off

The Python-side plumbing now emits and validates explicit JSON contracts before
real source decode is wired:

- casebuilder emits:
  - `prepared_case_manifest.json`
  - `analysis_state_stub.json`
  - `boundary_cache_stub.json`
  - `checkpoint_stub.json`
  - `product_plan.json`
- verification consumes those artifacts plus the current idealized dry summary
  JSON written by `gwm_idealized_driver`

That keeps the runtime core clean while still making the surrounding
source-to-product path concrete enough to integrate against.

## Near-Term Use

Near-term, the most useful reuse is:

1. `cfrust` for high-level HRRR / RRFS GRIB workflows
2. `ecrust` when low-level GRIB handle operations are needed
3. `wrf-rust` and `wrf-rust-plots` for output-side diagnostics and graphics

The new model should write outputs and metadata that are practical for those
downstream tools without distorting the internal runtime state model to mimic
WRF internals.
