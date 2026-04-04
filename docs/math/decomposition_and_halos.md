# Decomposition and Halos

## Continuous equations

There is no separate continuous PDE for decomposition. This module preserves
the meaning of the discrete dycore operators under a 2D horizontal tiling of
the global patch.

## Discrete update

- the global domain is partitioned into horizontal subdomains with full
  vertical columns on each rank/tile
- each subdomain owns halo-padded `Field3D` and `FaceField` storage
- scalar halo exchange is an ordered two-phase copy/transport:
  - pack west/east send faces
  - fill west/east ghost cells
  - pack south/north send faces
  - fill south/north ghost cells
- face halo exchange preserves the orientation-specific pack/unpack shapes for
  `mom_u`, `mom_v`, and `mom_w`
- X and Y owned-face interface synchronization is a separate step after halo
  exchange; Z-face synchronization is a no-op by contract
- the virtual-rank path is the reference implementation for pack order, absent
  neighbor behavior, and face-interface semantics
- the MPI milestone must preserve the same discrete contract while replacing
  only the transport backend inside `comm/*`

## Invariants / admissibility

- decomposition is a complete cover of the owned global domain with no overlap
  in owned-cell indices
- owned cells are never overwritten by unrelated halo copies
- halo width, face counts, and pack/unpack extents match the local tile
  geometry exactly
- periodic closure is respected when enabled
- scalar and face halo semantics are backend-invariant:
  - virtual-rank and MPI must produce the same ghost values for the same local
    layout
  - X/Y owned-face synchronization must produce the same post-sync interface
    state
- backend awareness is confined to `comm/*`; dycore, surface, and physics
  operators consume local haloed fields only

## Assumptions for stability / consistency

- all operator stages that need neighbor data perform halo exchange before
  reading ghost cells
- pack/unpack semantics use consistent face orientation and stride rules
- if MPI transport is active, Cartesian rank ordering matches the
  `VirtualRankLayout` row-major layout used by the oracle tests
- absent-neighbor semantics are explicit and backend-invariant; no backend may
  invent a different boundary fill rule

## Scope boundary

- topology discovery, halo plans, transport, and owned-face synchronization
  live in `comm/*`
- dycore code may request halo exchange or backend-neutral collectives, but it
  does not inspect ranks, communicators, or launcher state
- the MPI milestone does not widen `DryState`, change operator formulas, or
  introduce gather/step/scatter shortcuts

## Test mapping

Current local oracle coverage:

- `tests/unit/test_halo_exchange.cpp`
- `tests/unit/test_face_halo_exchange.cpp`
- `tests/unit/test_scalar_halo_buffers.cpp`
- `tests/unit/test_mpi_runtime.cpp`
- `tests/property/test_virtual_rank_equivalence.cpp`
- `tests/property/test_dry_virtual_rank_equivalence.cpp`
- `tests/property/test_dry_momentum_virtual_rank_equivalence.cpp`
- `tests/property/test_dry_fastmode_virtual_rank_equivalence.cpp`
- `tests/property/test_terrain_virtual_rank_equivalence.cpp`

Expected MPI milestone coverage:

Unit:

- `tests/unit/test_mpi_runtime.cpp`
- `tests/unit/test_mpi_cartesian_context.cpp`
- `tests/unit/test_mpi_scalar_halo_exchange.cpp`
- `tests/unit/test_mpi_face_halo_exchange.cpp`
- `tests/unit/test_mpi_face_interface_sync.cpp`

Property:

- `tests/property/test_mpi_halo_virtual_rank_equivalence.cpp`
- `tests/property/test_mpi_dry_virtual_rank_equivalence.cpp`
- `tests/property/test_mpi_dry_momentum_virtual_rank_equivalence.cpp`
- `tests/property/test_mpi_dry_fastmode_virtual_rank_equivalence.cpp`
- `tests/property/test_mpi_terrain_virtual_rank_equivalence.cpp`

Regression placeholders for the first real multi-rank dry gate:

- `tests/regression/test_mpi_hydrostatic_rest.cpp`
- `tests/regression/test_mpi_terrain_hydrostatic_rest.cpp`
- `tests/regression/test_mpi_density_current_evolution.cpp`
- `tests/regression/test_mpi_mountain_wave_evolution.cpp`
