# Milestones

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
- buoyancy-response regression for warm-bubble / density-current dry evolution
- horizontal-pressure-response regression for centered thermal perturbations
- dry diagnostics summary helpers
- native idealized dry driver plus casebuilder wrapper for idealized YAML cases
- next:
  - real MPI-backed halo exchange path once local/cluster MPI is available
  - terrain/hybrid-metric-aware dry operators on top of the now fuller dry core
  - meaningful mountain-wave dry regression with the fast path active

## stage2_moist

- land/surface state
- Monin-Obukhov surface layer
- TKE PBL
- warm-rain microphysics
- minimal radiation

## stage3_practical

- HRRR and RRFS adapters
- casebuilder target-area workflow
- severe-weather diagnostics/product set
- restart/rerun loop

## stage4_scaling

- real MPI halo exchange
- async output hardening
- 2-8 GPU tuning
- regression corpus for scaling and restart fidelity
