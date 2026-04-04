# MPI Halo Exchange

## Milestone scope

This note defines the bounded MPI transport milestone for the existing dry
core. The goal is to replace the current same-process virtual-rank transport
with a real rank-to-rank backend while preserving the established halo and
face-interface semantics.

The milestone is successful only if:

- MPI/backend awareness remains confined to `comm/*`
- the current virtual-rank path remains the oracle/reference implementation
- non-MPI builds remain behaviorally unchanged
- MPI-enabled runs match the existing local oracle within the documented
  tolerances

Out of scope:

- moist physics
- surface/land evolution
- ingest or boundary-cache work
- restart redesign
- map-factor or Coriolis terms
- changes to the dry operator formulas themselves

## Continuous equations

There is no separate continuous PDE for MPI transport. MPI halo exchange is a
backend realization of the existing discrete communication contract for
halo-padded fields on a 2D Cartesian decomposition.

## Discrete update

For a local subdomain with scalar field `q` and halo width `h`:

1. build or validate a Cartesian communication context whose rank ordering and
   periodic flags match the existing `VirtualRankLayout` contract
2. pack west/east scalar faces into explicit send buffers
3. exchange west/east buffers with direct neighbor ranks
4. unpack received west/east buffers into ghost cells
5. pack south/north scalar faces
6. exchange south/north buffers with direct neighbor ranks
7. unpack received south/north buffers into ghost cells

For face fields, the same ordered transport is applied using the
orientation-specific pack/unpack extents already defined by `FaceField`
geometry.

After halo exchange:

- X-face and Y-face owned interfaces are synchronized across rank boundaries
  using the same semantics as the local oracle path
- Z-face synchronization remains a no-op by contract

## Invariants / admissibility

- topology equivalence:
  - the active MPI Cartesian layout matches the row-major
    `VirtualRankLayout` ordering used by the existing oracle tests
- transport exactness:
  - scalar and face halos are filled with the same values the virtual-rank
    oracle would produce for the same decomposition
- owned-region preservation:
  - halo exchange never mutates unrelated owned cells
- absent-neighbor preservation:
  - missing-neighbor halo regions receive the same explicit fill behavior as
    the local oracle path
- interface consistency:
  - X/Y owned-face synchronization produces the same post-sync interface values
    as the oracle path
- backend confinement:
  - MPI symbols, communicators, rank logic, and launcher assumptions stay in
    `comm/*`

## Assumptions for stability / consistency

- communicator size equals `ranks_x * ranks_y`
- periodic flags and neighbor identities match the local decomposition
- the pack/unpack order remains X phase first, then Y phase
- any domain-wide reductions needed by the dry core are exposed through
  backend-neutral comm helpers rather than direct MPI use in dycore code
- test-only oracle comparisons may gather reference data, but production MPI
  transport may not use gather/step/scatter shortcuts

## Proof obligations

1. MPI Cartesian topology is equivalent to the existing local decomposition
   contract.
2. Scalar halo exchange is exact relative to the virtual-rank oracle.
3. Face halo exchange is exact relative to the virtual-rank oracle.
4. X/Y face-interface synchronization is exact relative to the oracle; Z
   remains unchanged.
5. Existing dry-core virtual-rank equivalence guarantees remain true when the
   backend is real MPI.
6. Non-MPI builds and MPI-enabled one-rank runs remain behaviorally unchanged.
7. No runtime layer outside `comm/*` becomes MPI-aware.

## Test mapping

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

## Acceptance shape

- local/non-MPI validation remains green with no tolerance loosening
- MPI-enabled one-rank runs reproduce current local behavior
- MPI-enabled multi-rank runs match the current oracle tolerances for scalar
  exchange, full dry stepping, momentum transport, fast modes, and terrain
  cases
- no production code path performs root gather, serial stepping, and scatter
- no file outside `comm/*` requires `mpi.h` or rank-aware logic
