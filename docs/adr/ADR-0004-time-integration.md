# ADR-0004: Time Integrator Family

## Status

Accepted.

## Decision

- SSPRK3 is the slow-term outer integrator.
- A local split-explicit fast-mode path handles fast acoustic/vertical
  processes.
- No global elliptic solve is included in the MVP.

## Rationale

This is a practical compromise between numerical discipline and GPU-first
execution without taking on a global pressure-solve dependency too early.

## Consequences

- The first dry core is organized around operator families and explicit stage
  storage.
- Fast-mode infrastructure is wired in early, even if the MVP implementation is
  still a scaffold.
