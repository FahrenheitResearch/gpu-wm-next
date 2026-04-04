# Prepared-Case Runtime

This note captures the JSON contracts that bridge source data into runtime
initialization and source-driven verification.

## Balanced Startup Conversion

Prepared-case dry-state initialization is not a raw source-field import.
Instead, ingest projects the source atmosphere onto the runtime dry manifold:

- compute source `theta` from source `air_pressure` and `air_temperature`
- integrate each column hydrostatically on `GridMetrics::z_center(...)` using
  the same four-iteration recurrence as the dry fast-mode reference builder
- write `rho_d` and `rho_theta_m` from that projected pressure column
- iteratively project the horizontal face momenta so the column-mean
  continuity residual is driven near roundoff before rebuilding `mom_w`
- build `mom_u` and `mom_v` from projected density plus face-interpolated
  source `u_wind` / `v_wind`
- initialize `mom_w` to a startup-consistent zero field with zero top/bottom
  faces; do not trust imported `w_wind` as prognostic startup truth
- rebuild tracer masses only after final projected `rho_d`

Prepared-case boundary snapshots use the same balanced conversion before the
open-boundary strips are applied.

The current boundary runtime path is intentionally strip-only:

- each rounded stage time caches one balanced boundary projection
- dry and tracer boundary updaters reuse that cached stage snapshot instead of
  rebuilding per-tracer whole-domain reference states
- application writes only the west/east/south/north strips that the open
  boundary path actually consumes

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

For the current warm-rain runtime milestone, the prepared-case bridge is
expected to preserve `specific_humidity` as the canonical runtime/product
name while allowing the runtime tracer set to grow to:

- `specific_humidity`
- `cloud_water_mixing_ratio`
- `rain_water_mixing_ratio`

For deeper source-driven columns, the populate step may:

- expand a legacy short pressure-level request into a deeper source-side
  pressure stack before vertical remap
- preserve the manifest `z_top` as the runtime target top instead of silently
  inflating it from source heights
- precondition supersaturated source humidity into initial
  `cloud_water_mixing_ratio` at the populate edge so the runtime does not take
  the entire condensation shock in the first timestep

The `water_vapor_mixing_ratio` alias remains acceptable only at the source
decode / ingest edge while the bridge is still being tightened.

When a companion `analysis_state.json` is present next to a source-driven
`plan_view.json`, the current product bridge may enrich the plan-view bundle
with moisture diagnostics derived from the populated analysis humidity field.

## Invariants / Admissibility

- populated analysis and boundary artifacts must match the declared grid shape
- field lengths must equal the corresponding 2-D or 3-D cell counts
- boundary offsets must be strictly increasing
- map manifests must reference existing image artifacts
- plan-view bundles must preserve row-major storage and exact field sizes
- moist enrichment fields, when present, must preserve the same plan-view slice
  shape and naming contract as dry fields
- prepared-case `summary.json` should remain finite and now includes:
  - dry totals/ranges
  - tracer totals/ranges
  - vapor/cloud/rain/condensed/total water masses
- startup balance diagnostics:
  - EOS residual against the balanced startup pressure column
  - hydrostatic residual against the balanced projected pressure field
  - fast-vertical residual from the perturbation state seen by the vertical
    fast-mode branch
  - mass-divergence residual from the startup face momenta
  - top/bottom `mom_w` norms
  - tracer closure
  - optional source-height vs metric-height mismatch

## Assumptions for Stability / Consistency

- source decode remains outside the runtime core
- the populated analysis/boundary artifacts are the runtime truth inputs
- the humidity alias is transitional and accepted only during the bridge period
- canonical runtime and product outputs should prefer `specific_humidity`
  even when the ingest bridge accepted an alias on input
- warm-rain summary sidecars are proof/report artifacts, not restart truth
- the source-run bundle should be stable enough to support first real maps
  before full moist physics lands
- populated metadata should record both the requested and actual source
  pressure-level stack used for remap
- balanced startup diagnostics are proof obligations for prepared-case ingest
  and must stay within the unit-test acceptance gates before longer real-data
  runs are trusted

## Test Mapping

- `tests/unit/test_runtime_case.cpp`
- `tests/unit/test_ingest_contracts.cpp`
- `tests/unit/test_prepared_case_init.cpp`
- `tests/unit/test_tracer_registry.cpp`
- `tests/unit/test_runtime_summary.cpp`
- `tests/regression/test_warm_rain_summary_closure.cpp`
- `tools/verify/run_verification.py`
- `tools/verify/test_source_run_bundle.py`
- `tools/verify/test_populate_prepared_case.py`
- `tools/verify/render_plan_view_maps.py`
