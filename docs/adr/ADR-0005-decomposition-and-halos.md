# ADR-0005: Decomposition and Halo Contracts

## Status

Accepted.

## Decision

- The model uses 2D horizontal decomposition from day one.
- Each rank owns complete vertical columns.
- Every field allocation includes halo cells.
- Single-GPU mode is a 1x1 decomposition, not a separate architecture.

## Rationale

This prevents later rewrites when moving from one GPU to several GPUs and keeps
physics column-local.

## Consequences

- A virtual-rank harness is required before real multi-rank deployment.
- Halo pack/unpack and exchange semantics are part of the core API.
- The implementation path should expose:
  - cartesian neighbor lookup
  - explicit face-buffer packing/unpacking
  - an MPI runtime boundary that degrades cleanly when MPI is unavailable
