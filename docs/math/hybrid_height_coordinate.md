# Hybrid Height Coordinate

## Continuous equations

The vertical coordinate is height-based near the top and terrain-following near
the surface:

`z(x, y, eta) = a(eta) * z_top + b(eta) * terrain_dyn(x, y)`

with `b(eta)` decreasing rapidly aloft.

## Discrete update

Layer interfaces and centers are precomputed from `a_k`, `b_k`, and smoothed
dynamic terrain. Metric terms are stored in `GridMetrics`.

## Invariants / admissibility

- monotone layer interfaces
- positive layer thickness
- exact reproduction of flat-terrain levels when terrain is zero

## Assumptions for stability / consistency

- terrain taper profile is smooth in `eta`
- dynamic terrain is smoothed enough to avoid pathological metric noise

## Test mapping

- future coordinate monotonicity/unit tests
- future flat-terrain equivalence tests
