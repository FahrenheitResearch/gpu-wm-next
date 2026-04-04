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

## Invariants / admissibility

- output schema version must be explicit
- each field must satisfy `len(values) == nx * ny`
- slice index must remain within `[0, nz - 1]`
- gathered plan-view fields must preserve constant-state values
- rendered map manifests must reference existing image artifacts

## Assumptions for stability / consistency

- current plan-view products are for the final dry state of a run
- current cell-centered velocity diagnostics are visualization-friendly derived
  products, not canonical prognostic state
- current rendering is idealized-case focused and does not yet imply real-data
  map parity

## Test mapping

- `tests/unit/test_plan_view_output.cpp`
- `tools/verify/run_verification.py`:
  - `plan_view_bundle`
  - `map_manifest`
- end-to-end smoke path through:
  - `tools/casebuilder/run_idealized_case.py`
  - `tools/verify/render_plan_view_maps.py`
