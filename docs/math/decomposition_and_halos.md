# Decomposition and Halos

## Continuous equations

There is no separate continuous equation. This module preserves the discrete
operator meaning under 2D horizontal decomposition.

## Discrete update

- the global domain is partitioned into horizontal tiles
- each tile owns a halo-padded local field
- neighbor-owned edge values are copied into ghost cells before stencil work
- x-face exchange and y-face exchange are treated as separate ordered phases
- explicit face buffers define the discrete communication contract

## Invariants / admissibility

- serial vs virtual-rank equivalence for supported operators
- owned cells are never overwritten by unrelated halo copies
- periodic closure is respected when enabled
- face-buffer sizes match halo width and local tile geometry exactly

## Assumptions for stability / consistency

- decomposition is a complete cover of the owned global domain
- pack/unpack semantics use consistent face orientation and stride rules

## Test mapping

- `tests/unit/test_halo_exchange.cpp`
- `tests/unit/test_face_halo_exchange.cpp`
- `tests/unit/test_scalar_halo_buffers.cpp`
- `tests/unit/test_mpi_runtime.cpp`
- `tests/property/test_virtual_rank_equivalence.cpp`
- `tests/property/test_dry_virtual_rank_equivalence.cpp`
- `tests/property/test_dry_momentum_virtual_rank_equivalence.cpp`
- `tests/property/test_dry_fastmode_virtual_rank_equivalence.cpp`
- `tests/property/test_terrain_virtual_rank_equivalence.cpp`
