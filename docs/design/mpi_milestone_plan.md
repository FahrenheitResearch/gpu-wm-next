# MPI Milestone Plan

This is a bounded implementation brief for the likely next runtime-core
milestone after `12aeedb`. It is not an architecture reopen.

## Goal

Replace the current MPI runtime stub with a real MPI-backed halo exchange path
for the existing dry core, while preserving:

- the current single-state GPU-first architecture
- the existing virtual-rank oracle semantics
- the rule that only `comm/*` is MPI/backend aware

## In scope

- real MPI Cartesian topology / neighbor mapping
- scalar halo exchange for `Field3D`
- face-field halo exchange for `FaceField`
- shared X/Y face-interface synchronization across rank boundaries
- MPI-vs-virtual-rank equivalence tests for:
  - scalar halo exchange
  - full dry SSPRK3 stepping
  - dry momentum transport
  - fast modes
  - terrain-aware dry cases

## Out of scope

- moist physics
- surface/land flux evolution
- ingest/boundary-cache work
- restart redesign
- Coriolis, damping, diffusion, map-factor work
- terrain-operator changes
- MPI logic outside `comm/*`

## Acceptance shape

The milestone is only complete if:

- non-MPI builds remain green
- MPI-enabled 1-rank runs match existing local behavior
- MPI-enabled multi-rank runs match the current virtual-rank oracle within the
  established property/regression tolerances
- no operator code in `dycore/*`, `surface/*`, or `physics/*` gains MPI/backend
  awareness

## Likely file touch list

- `include/gwm/comm/halo_exchange.hpp`
- `src/comm/halo_exchange.cpp`
- `include/gwm/comm/mpi_runtime.hpp`
- `src/comm/mpi_runtime.cpp`
- `docs/math/decomposition_and_halos.md`
- `CMakeLists.txt`

Likely new files:

- `docs/math/mpi_halo_exchange.md`
- `cmake/GwmAddMpiTest.cmake`
- `tests/unit/test_mpi_cartesian_topology.cpp`
- `tests/unit/test_mpi_scalar_halo_exchange.cpp`
- `tests/unit/test_mpi_face_halo_exchange.cpp`
- `tests/unit/test_mpi_halo_mismatch_rejects.cpp`
- `tests/property/test_mpi_virtual_rank_equivalence.cpp`
- `tests/property/test_mpi_dry_virtual_rank_equivalence.cpp`
- `tests/property/test_mpi_dry_momentum_virtual_rank_equivalence.cpp`
- `tests/property/test_mpi_dry_fastmode_virtual_rank_equivalence.cpp`
- `tests/property/test_mpi_terrain_virtual_rank_equivalence.cpp`
- `tests/regression/test_mpi_hydrostatic_rest.cpp`
- `tests/regression/test_mpi_terrain_hydrostatic_rest.cpp`
- `tests/regression/test_mpi_density_current_evolution.cpp`
- `tests/regression/test_mpi_mountain_wave_evolution.cpp`

## Guardrails

- Do not widen `DryState` just to make MPI easier.
- Do not move neighbor lookup or pack/unpack semantics out of `comm/*`.
- Do not loosen existing dry-core tolerances to make MPI look acceptable.
- Do not add host-side gather/step/scatter shortcuts.
- Do not mix this milestone with new physics or surface work.
