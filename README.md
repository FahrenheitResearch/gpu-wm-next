# gpu-wm-next

Fresh-start GPU-first mesoscale weather model scaffold.

This repository is the successor architecture build wave, not a refactor of
`gpu-wm`. The old model remains a reference/oracle for:

- ideas that worked
- failures that should not be repeated
- verification workflows
- visual/diagnostic expectations

The implementation priorities for this repository are:

1. lock the core architecture early
2. keep the runtime GPU-native
3. build multi-GPU-safe data structures and contracts from day one
4. ship proof-obligation docs and targeted tests alongside code

Start with:

- [AGENTS.md](AGENTS.md)
- [docs/design/milestones.md](docs/design/milestones.md)
- [docs/adr](docs/adr)
- [docs/math](docs/math)
- [docs/design/external_toolchain_integration.md](docs/design/external_toolchain_integration.md)

Local Windows/CUDA build:

- `powershell -ExecutionPolicy Bypass -File scripts/build_local_windows.ps1`

Run an idealized dry case locally:

- `python tools/casebuilder/run_idealized_case.py cases/idealized/warm_bubble.yaml --steps 4 --dt 2.0 --fast-substeps 3`

Current implementation status:

- Stage 0 scaffold is in place and documented.
- Stage 1 dry now includes:
  - rectangle-domain builder
  - topography-aware idealized domain builder support
  - halo-aware dry conserved-state bundle
  - SSPRK3 dry-state scaffold
  - local split-explicit fast-mode acoustic subcycle
  - EOS-consistent dry pressure reconstruction from `rho_d` and
    `rho_d*theta_m`
  - conservative first-order slow momentum-flux transport for `mom_u`,
    `mom_v`, and `mom_w`
  - face-halo exchange for C-grid momentum fields
  - horizontal pressure-gradient momentum coupling for `mom_u` and `mom_v`
  - perturbation-form vertical fast-mode treatment plus buoyancy-driven
    vertical-momentum slow tendency in `mom_w`
  - hydrostatic-rest and dry-bundle regression binaries
  - buoyancy-response and horizontal-pressure-response regressions for
    warm-bubble and density-current style forcing
  - density-current evolution regression on the pressure-coupled fast/slow dry
    core
  - acoustic-pulse regression and fast-mode virtual-rank equivalence tests
  - terrain/hybrid-height metric-aware dry operators with flat-limit coverage
  - terrain hydrostatic-rest and mountain-wave evolution regressions
  - verified additive composition of slow momentum sources
  - fixed mountain-wave conserved-momentum initialization
  - idealized dry state initializers for warm-bubble, density-current, and
    mountain-wave background cases
  - dry diagnostics summary helpers for mass/theta/momentum/extrema
  - native idealized driver executable (`gwm_idealized_driver`)
  - Python casebuilder wrapper for idealized YAML case files
  - prepared-case/source-driven runtime bridge for populated
    `analysis_state.json` and `boundary_cache.json` artifacts
  - prepared-case runtime driver executable (`gwm_prepared_case_driver`)
  - Python source-run wrapper for prepared-case manifests
  - first actual-data smoke path from HRRR into rendered plan-view map products
  - source-run bundle verification for populated analysis, boundary, summary,
    plan-view, and map-manifest JSON artifacts
- external-tooling boundary is now codified for:
  - `cfrust`
  - `ecrust`
  - `wrf-rust`
  - `wrf-rust-plots`

Windows note:

- Full local validation now runs through `ctest` directly under
  `scripts/build_local_windows.ps1`.
- `scripts/run_gate_b_windows.cmd` remains available as a targeted Stage 1 dry
  regression wrapper.

Current checkpoint:

- The local non-MPI suite passes with the terrain-aware fast/slow dry core and
  prepared-case runtime path in place.
- Real MPI-backed halo exchange is implemented source-side and covered by the
  comm-layer proof/docs; execution of the MPI tests is still pending a machine
  with `mpiexec`.
- The repo now has two working map paths:
  - idealized dry cases through `gwm_idealized_driver`
  - populated prepared/source-driven runs through `gwm_prepared_case_driver`

Source-driven smoke example:

- `python tools/casebuilder/prepare_case.py --source hrrr --domain-name actual_data_smoke --cycle-time-utc 2026-04-04T00:00:00Z --forecast-hours 1 --nx 24 --ny 24 --center-lat 39.0 --center-lon -97.0 --pressure-levels-hpa 1000,850,700`
- `python tools/casebuilder/run_prepared_case.py --prepared-case cases/prepared/actual_data_smoke_hrrr/prepared_case_manifest.json --populate --steps 1 --dt 2.0 --fast-substeps 2 --plan-view-level 1 --render-maps`
- `python tools/verify/run_verification.py --input cases/prepared/actual_data_smoke_hrrr --kind source_run_bundle`
