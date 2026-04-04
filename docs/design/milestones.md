# Milestones

## execution_mode

- default mode from this checkpoint forward: direct implementation
- do not stop for routine GPT Pro checkpoints
- use GPT Pro only for major architecture forks or expensive-to-reverse
  milestone cuts
- current checkpoint for that policy: `b2dfdad`

## stage0

- buildable CMake/CUDA scaffold
- ADR set for the frozen architecture choices
- math notes for state, coordinates, dycore, decomposition, screen obsops,
  surface fluxes
- field/state/decomposition abstractions with tests

## stage1_dry

- dry transport kernel path
- SSPRK3 slow-step skeleton
- fast-mode integrator interface
- halo-aware virtual-rank harness
- constant-state and conservation regression cases
- rectangle-domain builder
- dry conserved-state bundle (`rho_d`, `rho_d*theta_m`, face momentum)
- hydrostatic-rest placeholder regression
- EOS-consistent dry pressure reconstruction from conservative state
- first-order conservative slow momentum-flux transport on the flat-grid dry
  core
- face-halo exchange for C-grid momentum fields
- horizontal pressure-gradient momentum coupling for `mom_u` and `mom_v`
- local split-explicit acoustic subcycle with perturbation-form vertical update
- buoyancy-driven vertical-momentum response for thermal perturbations
- cartesian-neighbor / face-buffer comm contract
- MPI runtime boundary stub with clean no-MPI behavior
- idealized dry case library for warm-bubble / density-current / mountain-wave
  background states
- topography-aware idealized domain builder support
- terrain/hybrid-height metric-aware dry operators threaded through the fast and
  slow dry core
- terrain-aware hydrostatic-rest and mountain-wave dry regressions with the fast
  path active
- fixed slow-source composition semantics across momentum-flux,
  pressure-gradient, and buoyancy tendencies
- fixed mountain-wave background-flow initialization into conserved face
  momentum
- buoyancy-response regression for warm-bubble / density-current dry evolution
- horizontal-pressure-response regression for centered thermal perturbations
- dry diagnostics summary helpers
- native idealized dry driver plus casebuilder wrapper for idealized YAML cases
- next:
  - real MPI-backed halo exchange path against the existing virtual-rank oracle
  - keep `comm/*` as the only MPI/backend-aware runtime layer
  - use `docs/math/mpi_halo_exchange.md` as the bounded proof/scope note for
    the transport milestone
  - keep surface/land spine and ingest/boundary-cache work gated behind the
    current milestone decision

## stage2_moist

- land/surface state
- minimal surface runtime spine:
  - `SurfaceState`
  - `SurfaceStaticProperties`
  - shared neutral surface-layer closure
  - tile aggregation helper for screen diagnostics / future fluxes
- Monin-Obukhov surface layer
- TKE PBL
- warm-rain microphysics
- minimal radiation

## stage3_practical

- HRRR and RRFS adapters
- canonical ingest/runtime schema for analysis state and time-indexed boundary
  cache metadata
- casebuilder target-area workflow
- severe-weather diagnostics/product set
- restart/rerun loop

## stage4_scaling

- multi-rank dry regressions and scaling hardening
- async output hardening
- 2-8 GPU tuning
- regression corpus for scaling and restart fidelity
