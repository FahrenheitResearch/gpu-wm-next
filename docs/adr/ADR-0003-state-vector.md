# ADR-0003: Canonical Prognostic State

## Status

Accepted.

## Decision

The atmospheric prognostic state is:

- `rho_d`
- `rho_d * u`
- `rho_d * v`
- `rho_d * w`
- `rho_d * theta_m`
- tracer masses managed by a registry

The surface state is separate and tile-capable.

## Rationale

A conservative mass-based state fits finite-volume transport, tracer management,
and later microphysics coupling while staying GPU-friendly.

## Consequences

- All operators consume/produce the canonical state or explicit tendencies.
- Screen-level diagnostics are observation operators, not alternate truth
  fields.
