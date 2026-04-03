# ADR-0008: External Tooling Boundary

## Status

Accepted.

## Decision

- The runtime core remains one C++20 + CUDA + MPI codebase.
- External Rust tooling is a first-class part of the development stack, but it
  stays outside the timestepper / physics runtime.
- Preferred external tooling roles are:
  - `cfrust`: high-level GRIB ingest
  - `ecrust`: low-level GRIB/ecCodes-style escape hatch
  - `wrf-rust`: diagnostics / verification
  - `wrf-rust-plots`: plotting / products
  - Rust MetPy replacement repo: generic meteorological math
  - `rusbie` / `rustbie`: acquisition workflow

## Rationale

This keeps the solver architecture clean while still exploiting the user's
existing weather-tooling ecosystem. It accelerates case preparation,
verification, and product generation without creating a split-runtime core.

## Consequences

- The ingest layer should expose clean canonical-IR boundaries where external
  tools can plug in.
- Tooling scripts should know how to discover and use those repos.
- No solver/physics code should assume those repos are available at runtime for
  timestep execution.
