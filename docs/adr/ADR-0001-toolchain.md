# ADR-0001: Toolchain and Runtime Language

## Status

Accepted.

## Decision

- Runtime core uses C++20 + CUDA + CMake.
- MPI is the communication interface for multi-GPU execution.
- Python is restricted to tooling, verification, plotting, and offline case
  helpers.
- No portability layer is used in the runtime core.

## Rationale

The project is explicitly CUDA-first. The runtime should optimize for the local
RTX 5090 and later rented A100 nodes rather than preserving CPU portability.

## Consequences

- Runtime kernels and memory layout are written directly against CUDA.
- MPI hooks exist from day one, even before cluster-scale tests.
- Rust/Python side projects may support ingest, plotting, and verification, but
  not own the timestepper.
