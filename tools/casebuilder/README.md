# Casebuilder

This directory will host the target-area domain builder, source-ingest wrappers,
and boundary-cache generation tools.

Stage 0 keeps only the contract:

- source-specific decode belongs here
- canonical IR generation belongs here
- runtime core does not parse raw source formats directly

Preferred surrounding tooling:

- `rusbie` / `rustbie` for acquisition workflows
- `cfrust` for primary GRIB decode
- `ecrust` for low-level GRIB fallback

Current helper entry points:

- `check_external_toolchain.py`
- `prepare_case.py`
- `run_idealized_case.py`
- `toolchain_manifest.example.toml`
