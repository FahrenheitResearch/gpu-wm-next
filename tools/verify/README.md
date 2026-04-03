# Verification

This directory will host the successor-model verification harness.

Core principles carried forward from the prototype:

- fixed sample sets
- candidate vs control comparisons
- terrain-aware evaluation
- screen-level diagnostics vs optional proxy comparisons
- fast falsification over aesthetic branch drift

Preferred surrounding tooling:

- `wrf-rust` for diagnostics / verification
- `wrf-rust-plots` for product generation
- Rust MetPy replacement repo for generic thermo / severe-weather math

Current helper entry point:

- `run_verification.py`
