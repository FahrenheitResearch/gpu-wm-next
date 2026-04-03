# Project Mission

Build a new GPU-first mesoscale weather model from scratch. Do not patch or
extend `gpu-wm` except as a reference for ideas, tests, or comparisons.

# Non-Negotiable Architecture

- Core language/toolchain: C++20 + CUDA + CMake + MPI.
- Python is allowed only for tooling, verification scripts, plotting, and
  offline helpers.
- Rust tooling is explicitly encouraged around the runtime core when it
  accelerates ingest, diagnostics, plotting, or workflow automation.
- No Kokkos/SYCL/portability abstraction in the runtime core.
- No CPU shadow state for the model core.
- No hidden single-GPU assumptions.
- No WRF/MPAS code import or architecture copying.
- No surface humidity endpoint blending, terrain-gated moisture hacks, or
  LML-proxy diagnostics masquerading as truth.

# External Tooling Boundary

- Keep the timestepper, dynamics, physics, surface model, decomposition, and
  restart state in one C++20/CUDA/MPI runtime.
- Use surrounding Rust tooling aggressively for:
  - model-data acquisition
  - GRIB ingest helpers
  - diagnostics / verification
  - plotting / product generation
  - case preparation workflows
- Preferred surrounding repos from prior evaluation:
  - `cfrust`
  - `ecrust`
  - `wrf-rust`
  - `wrf-rust-plots`
  - Rust MetPy replacement repo when pinned
  - `rusbie` / `rustbie` when pinned
- Do not make the solver core depend on those repos for timestep execution.

# Numerical Choices Already Fixed

- Structured projected regional grid; Lambert conformal first.
- Hybrid terrain-following height coordinate with terrain taper aloft.
- Finite-volume compressible nonhydrostatic core on C-grid staggering.
- State variables:
  - `rho_d`
  - `rho_d*u`
  - `rho_d*v`
  - `rho_d*w`
  - `rho_d*theta_m`
  - tracer masses managed by a registry
- Time stepping:
  - SSPRK3 for slow terms
  - local split-explicit fast-mode path
  - no global elliptic solve in the MVP
- Multi-GPU architecture:
  - 2D horizontal decomposition
  - halo/ghost-cell API from day one
  - virtual-rank decomposition harness before actual cluster runs

# Surface / Physics Choices Already Fixed

- Tile-capable land model from day one; MVP can use `ntile=1`.
- 4 soil layers, canopy water bucket, simple snow branch, water branch.
- Monin-Obukhov surface layer and screen observation operators.
- PBL: one local TKE scheme first.
- Microphysics: warm-rain first, then one single-moment mixed-phase scheme.
- Radiation: simple broadband/two-stream first.

# Proof Obligations

Every new module/operator must include:

1. `docs/math/<module>.md`
2. governing equations
3. discrete update equations
4. stated invariants / admissibility conditions
5. assumptions for stability / consistency
6. exact tests that verify the stated properties

Minimum examples:

- dycore: conservation, constant-state preservation, hydrostatic-rest test
- tracers: positivity
- decomposition: serial vs virtual-rank equivalence
- screen obsops: analytic Monin-Obukhov limit tests
- microphysics: closed-box mass conservation
- radiation: clear-sky and zero-optical-depth limits

# Working Mode

- Main Codex thread acts as architect/integrator, not a monolithic coder.
- Spawn subagents for parallelizable work.
- Main thread resolves conflicts, integrates, builds, and runs final tests.
- Do not ask for confirmation on routine engineering decisions.
- Stop only for blockers that change equations, state vector, decomposition
  contract, or external dependency strategy.

# Deliverable Style

- Prefer small, buildable increments.
- Run tests after every meaningful patch set.
- Keep design docs in sync with code.
- When creating TODOs, group them by milestone:
  - `stage0`
  - `stage1_dry`
  - `stage2_moist`
  - `stage3_practical`
  - `stage4_scaling`
