from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare a future target-area case using external ingest tooling."
    )
    parser.add_argument("--source", required=True, choices=["hrrr", "rrfs", "rap", "gfs", "era5", "ecmwf"])
    parser.add_argument("--domain-name", required=True)
    parser.add_argument("--manifest", default="tools/casebuilder/toolchain_manifest.example.toml")
    args = parser.parse_args()

    manifest = Path(args.manifest)
    print(f"Casebuilder skeleton only. Source={args.source} domain={args.domain_name}")
    print(f"Toolchain manifest: {manifest.resolve()}")
    print("Expected future path:")
    print("- acquisition via rusbie/rustbie when available")
    print("- high-level GRIB decode via cfrust")
    print("- low-level escape hatch via ecrust")
    print("- canonical IR handoff into gpu-wm-next ingest layer")


if __name__ == "__main__":
    main()
