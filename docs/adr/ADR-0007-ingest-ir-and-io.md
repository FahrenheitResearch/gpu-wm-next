# ADR-0007: Canonical Ingest IR and I/O Strategy

## Status

Accepted.

## Decision

- All source datasets are converted into one canonical intermediate
  representation before they become model state.
- Source-specific adapters are thin and only decode/map source peculiarities.
- External Rust tooling is allowed and preferred around ingest/verification
  boundaries, but it does not enter the timestepper/runtime core.
- Output is split into:
  - native restart/checkpoint data
  - async NetCDF exports for analysis workflows

## Rationale

This keeps the model architecture independent of HRRR/RRFS/GFS/etc. quirks and
avoids forcing internal layout to mimic WRF files.

## Consequences

- Initialization and boundary generation are cleanly separable from the
  timestepper.
- Restart/checkpoint format can optimize for model fidelity instead of external
  compatibility.
