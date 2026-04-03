# ADR-0006: Surface Architecture and Screen Operators

## Status

Accepted.

## Decision

- Surface thermodynamics use a real energy-balance-driven surface state, not a
  minimal slab as the long-term architecture.
- Surface moisture must emerge from land/surface state and conductance, not
  heuristic endpoint blending.
- T2, Q2, RH2, U10, and V10 are first-class screen observation operators.

## Rationale

The previous prototype showed that honest screen operators were a real win and
that terrain-textured moisture/gate logic was not the right long-term path.

## Consequences

- Surface and obs-operator architecture are first-class design surfaces in the
  new model.
- Proxy model-level diagnostics may exist for verification only and must be
  labeled as such.
