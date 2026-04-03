# Hybrid Height Coordinate

## Continuous equations

The vertical coordinate is height-based near the top and terrain-following near
the surface:

`z(x, y, eta) = a(eta) * z_top + b(eta) * terrain_dyn(x, y)`

with `b(eta)` decreasing rapidly aloft.

## Discrete update

Layer interfaces and centers are precomputed from `a_k`, `b_k`, and smoothed
dynamic terrain. `GridMetrics` now stores terrain-aware interface/center
heights, inverse cell/face thicknesses, and terrain slopes, while preserving the
exact flat limit when `terrain_dyn = 0`.

## Invariants / admissibility

- monotone layer interfaces
- positive layer thickness
- exact reproduction of flat-terrain levels when terrain is zero

## Assumptions for stability / consistency

- terrain taper profile is smooth in `eta`
- dynamic terrain is smoothed enough to avoid pathological metric noise

## Test mapping

- `test_hybrid_height_metrics`
- `test_terrain_hydrostatic_rest`
- `test_terrain_virtual_rank_equivalence`
- `test_mountain_wave_evolution`
