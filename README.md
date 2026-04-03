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
  - local split-explicit fast-mode interface
  - buoyancy-driven vertical-momentum slow tendency for thermal perturbation
    response
  - hydrostatic-rest and dry-bundle regression binaries
  - buoyancy-response regression for warm-bubble and density-current forcing
  - idealized dry state initializers for warm-bubble, density-current, and
    mountain-wave background cases
  - dry diagnostics summary helpers for mass/theta/vertical-momentum/extrema
  - native idealized driver executable (`gwm_idealized_driver`)
  - Python casebuilder wrapper for idealized YAML case files
- external-tooling boundary is now codified for:
  - `cfrust`
  - `ecrust`
  - `wrf-rust`
  - `wrf-rust-plots`

Windows note:

- The new dry-state regression binaries currently run cleanly through the local
  wrapper scripts, and the Windows test wrapper forces synchronous CUDA launch
  so CTest sees deterministic exits. `scripts/build_local_windows.ps1` is the
  supported local validation path for now.
