# Prepared-Case Runtime

This note captures the JSON contracts that bridge source data into runtime
initialization and source-driven verification.

## Analysis-State Contract

The populated analysis state is a canonical runtime input with:

- `schema_version = "gwm-next-analysis-state/v1"`
- positive `nx`, `ny`, `nz`, `dx`, `dy`, `z_top`
- populated `atmosphere`, `surface`, and `static_surface` bundles
- populated metadata, including a `status` of `populated`
- a `target_window` block for the source-run workflow

Required atmosphere fields:

- `u_wind`
- `v_wind`
- `w_wind`
- `air_temperature`
- `air_pressure`
- `geopotential_height`
- a humidity field under either `specific_humidity` or the transitional
  `water_vapor_mixing_ratio` alias

Required surface fields:

- `surface_pressure`
- `air_temperature_2m`
- `specific_humidity_2m`
- `u_wind_10m`
- `v_wind_10m`
- `skin_temperature`

Required static-surface fields:

- `terrain_height`
- `land_mask`
- `land_use_index`

## Boundary-Cache Contract

The populated boundary cache is a time-indexed runtime input with:

- `schema_version = "gwm-next-boundary-cache/v1"`
- positive `boundary_interval_seconds`
- positive grid dimensions matching the analysis state
- at least two monotonically increasing snapshots
- snapshot-level atmosphere and surface bundles with the same required fields
  as the analysis state

The boundary cache is valid when each snapshot is individually populated and
the offsets are strictly increasing.

## Source-Run Bundle

The source-run bundle is a directory-level contract that ties together:

- `prepared_case_manifest.json`
- `analysis_state.json`
- `boundary_cache.json`
- `summary.json`
- `plan_view.json`
- `map_manifest.json`

This bundle is not a restart format. It is the source-to-runtime-to-product
bridge used for verification and map generation.

## Invariants / Admissibility

- populated analysis and boundary artifacts must match the declared grid shape
- field lengths must equal the corresponding 2-D or 3-D cell counts
- boundary offsets must be strictly increasing
- map manifests must reference existing image artifacts
- plan-view bundles must preserve row-major storage and exact field sizes

## Assumptions for Stability / Consistency

- source decode remains outside the runtime core
- the populated analysis/boundary artifacts are the runtime truth inputs
- the humidity alias is transitional and accepted only during the bridge period
- the source-run bundle should be stable enough to support first real maps
  before full moist physics lands

## Test Mapping

- `tests/unit/test_runtime_case.cpp`
- `tests/unit/test_ingest_contracts.cpp`
- `tools/verify/run_verification.py`
- `tools/verify/test_source_run_bundle.py`
- `tools/verify/render_plan_view_maps.py`
