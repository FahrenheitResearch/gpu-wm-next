# Terrain Metric Operators

## Continuous equations

The dry core uses a hybrid terrain-following height coordinate

`z(x, y, eta) = z_flat(eta) + b(eta) * h(x, y)`

with `h(x, y)` from `terrain_dyn` and `b(eta)` tapering to zero aloft. On this
coordinate, the horizontal pressure-gradient operator uses the chain-rule form

`∂p/∂x|z = ∂p/∂x|eta - (∂z/∂x|eta) * ∂p/∂z`

and likewise in `y`.

## Discrete update

- `GridMetrics` precomputes terrain-aware interface/center heights and inverse
  vertical spacings.
- Continuity, tracer/theta transport, fast acoustic density updates, and slow
  momentum transport consume local `inv_dz` fields rather than `dz_nominal`.
- Horizontal pressure-gradient operators consume terrain-aware center heights to
  apply a slope correction while retaining the exact flat-limit path when
  `terrain_dyn = 0`.

## Invariants / admissibility

- exact flat-terrain reduction to the prior Cartesian dry core
- monotone interface heights
- positive layer thickness everywhere
- no terrain-only source terms outside the metrics-driven operators

## Assumptions for stability / consistency

- terrain taper is smooth and bounded
- dynamic terrain remains smooth enough that discrete slope terms do not inject
  grid-scale noise
- shared C-grid interface faces are synchronized across virtual ranks before
  divergence updates

## Test mapping

- `test_hybrid_height_metrics`
- `test_terrain_virtual_rank_equivalence`
- `test_terrain_hydrostatic_rest`
- `test_mountain_wave_evolution`
- existing flat-limit regressions (`test_hydrostatic_rest`,
  `test_dry_fast_modes`, `test_dry_momentum_virtual_rank_equivalence`)
