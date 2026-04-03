# ADR-0002: Horizontal Grid and Vertical Coordinate

## Status

Accepted.

## Decision

- Use a structured, logically rectangular projected regional grid.
- Lambert conformal is the first supported projection.
- Use a hybrid terrain-following height coordinate with aggressive terrain taper
  aloft.

## Rationale

This matches the severe-weather target-area mission while keeping halo exchange,
column physics, and ingest straightforward.

## Consequences

- Domain decomposition is 2D horizontal only.
- Terrain enters both the dycore metrics and representativeness metadata, but
  those roles stay separate.
