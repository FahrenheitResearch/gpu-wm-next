from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    tomllib = None


REPO_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class FieldBinding:
    source_name: str
    canonical_name: str


@dataclass(frozen=True)
class SourceContract:
    key: str
    display_name: str
    boundary_interval_seconds: int
    first_class: bool
    atmosphere: tuple[FieldBinding, ...]
    surface: tuple[FieldBinding, ...]
    static_surface: tuple[FieldBinding, ...]


SOURCE_CONTRACTS: dict[str, SourceContract] = {
    "hrrr": SourceContract(
        key="hrrr",
        display_name="HRRR",
        boundary_interval_seconds=3600,
        first_class=True,
        atmosphere=(
            FieldBinding("UGRD", "u_wind"),
            FieldBinding("VGRD", "v_wind"),
            FieldBinding("VVEL", "w_wind"),
            FieldBinding("T", "air_temperature"),
            FieldBinding("SPFH", "specific_humidity"),
            FieldBinding("PRES", "air_pressure"),
            FieldBinding("HGT", "geopotential_height"),
        ),
        surface=(
            FieldBinding("PSFC", "surface_pressure"),
            FieldBinding("T2", "air_temperature_2m"),
            FieldBinding("Q2", "specific_humidity_2m"),
            FieldBinding("U10", "u_wind_10m"),
            FieldBinding("V10", "v_wind_10m"),
            FieldBinding("TSK", "skin_temperature"),
        ),
        static_surface=(
            FieldBinding("HGT", "terrain_height"),
            FieldBinding("LANDMASK", "land_mask"),
            FieldBinding("LU_INDEX", "land_use_index"),
        ),
    ),
    "rrfs": SourceContract(
        key="rrfs",
        display_name="RRFS",
        boundary_interval_seconds=3600,
        first_class=True,
        atmosphere=(
            FieldBinding("UGRD", "u_wind"),
            FieldBinding("VGRD", "v_wind"),
            FieldBinding("DZDT", "w_wind"),
            FieldBinding("TMP", "air_temperature"),
            FieldBinding("SPFH", "specific_humidity"),
            FieldBinding("PRES", "air_pressure"),
            FieldBinding("HGT", "geopotential_height"),
        ),
        surface=(
            FieldBinding("PRES", "surface_pressure"),
            FieldBinding("TMP_2M", "air_temperature_2m"),
            FieldBinding("SPFH_2M", "specific_humidity_2m"),
            FieldBinding("UGRD_10M", "u_wind_10m"),
            FieldBinding("VGRD_10M", "v_wind_10m"),
            FieldBinding("TMP_SFC", "skin_temperature"),
        ),
        static_surface=(
            FieldBinding("HGT_SFC", "terrain_height"),
            FieldBinding("LANDMASK", "land_mask"),
            FieldBinding("LANDUSE", "land_use_index"),
        ),
    ),
}

OPTIONAL_TOOL_ROLES = {
    "cfrust": "high-level GRIB ingest",
    "ecrust": "low-level GRIB/ecCodes-style access",
    "wrf_rust": "diagnostics and verification",
    "wrf_rust_plots": "plot rendering",
    "metrust": "thermodynamic/severe-weather math",
    "rusbie": "source-data acquisition workflow",
}


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def parse_utc(text: str) -> datetime:
    normalized = text.strip()
    if normalized.endswith("Z"):
        normalized = normalized[:-1] + "+00:00"
    return datetime.fromisoformat(normalized).astimezone(UTC)


def format_utc(value: datetime) -> str:
    return value.astimezone(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def default_cycle_time() -> str:
    # Use a conservative lag so the default prepared case targets a cycle that
    # is more likely to have finished indexing on remote object storage.
    now = (datetime.now(UTC) - timedelta(hours=2)).replace(
        minute=0, second=0, microsecond=0
    )
    return format_utc(now)


def parse_toolchain_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"paths": {}, "policy": {}, "exists": False}
    if tomllib is None:
        raise RuntimeError("Python tomllib is unavailable; cannot parse toolchain TOML")
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    data["exists"] = True
    return data


def canonical_names(bindings: tuple[FieldBinding, ...]) -> list[str]:
    return [binding.canonical_name for binding in bindings]


def binding_rows(bindings: tuple[FieldBinding, ...]) -> list[dict[str, str]]:
    return [
        {"source_name": binding.source_name, "canonical_name": binding.canonical_name}
        for binding in bindings
    ]


def build_offsets(interval_seconds: int, forecast_hours: int) -> list[int]:
    if forecast_hours < 0:
        raise ValueError("forecast-hours must be non-negative")
    forecast_seconds = forecast_hours * 3600
    offsets = list(range(0, forecast_seconds + interval_seconds, interval_seconds))
    if offsets[-1] != forecast_seconds:
        offsets.append(forecast_seconds)
    seen: set[int] = set()
    ordered: list[int] = []
    for offset in offsets:
        if offset not in seen:
            seen.add(offset)
            ordered.append(offset)
    return ordered


def parse_pressure_levels(text: str) -> list[int]:
    parts = [part.strip() for part in text.split(",") if part.strip()]
    if not parts:
        raise ValueError("pressure-levels-hpa must not be empty")
    levels = [int(part) for part in parts]
    if any(level <= 0 for level in levels):
        raise ValueError("pressure levels must be positive hPa values")
    return levels


def make_tool_records(tool_manifest: dict[str, Any]) -> list[dict[str, Any]]:
    configured_paths = tool_manifest.get("paths", {})
    tools: list[dict[str, Any]] = []
    for name, role in OPTIONAL_TOOL_ROLES.items():
        raw_path = configured_paths.get(name)
        resolved = Path(raw_path).expanduser().resolve() if raw_path else None
        tools.append(
            {
                "name": name,
                "role": role,
                "configured_path": str(resolved) if resolved else "",
                "exists": bool(resolved and resolved.exists()),
            }
        )
    return tools


def product_plan_stub(domain_name: str, source: SourceContract) -> dict[str, Any]:
    return {
        "schema_version": "gwm-next-product-plan/v1",
        "domain_name": domain_name,
        "source": source.display_name,
        "verification_modes": ["screen", "aloft", "products"],
        "current_supported_inputs": [
            {
                "kind": "prepared_case_manifest",
                "note": "Cold-start contract and boundary-cache planning",
            },
            {
                "kind": "idealized_summary_json",
                "note": "Current dry-core summary/product bridge from gwm_idealized_driver",
            },
        ],
        "planned_products": {
            "screen": [
                "surface_pressure",
                "skin_temperature",
                "air_temperature_2m",
                "specific_humidity_2m",
                "u_wind_10m",
                "v_wind_10m",
            ],
            "aloft": canonical_names(source.atmosphere),
            "products": [
                "verification_report_json",
                "plot_manifest_json",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare a concrete source-side case manifest and stub artifacts."
    )
    parser.add_argument(
        "--source",
        required=True,
        choices=["hrrr", "rrfs", "rap", "gfs", "era5", "ecmwf"],
    )
    parser.add_argument("--domain-name", required=True)
    parser.add_argument("--manifest", default="tools/casebuilder/toolchain_manifest.example.toml")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--cycle-time-utc", default=default_cycle_time())
    parser.add_argument("--forecast-hours", type=int, default=6)
    parser.add_argument("--nx", type=int, default=200)
    parser.add_argument("--ny", type=int, default=200)
    parser.add_argument("--nz", type=int, default=60)
    parser.add_argument("--dx", type=float, default=3000.0)
    parser.add_argument("--dy", type=float, default=3000.0)
    parser.add_argument("--z-top", type=float, default=20000.0)
    parser.add_argument("--center-lat", type=float, default=None)
    parser.add_argument("--center-lon", type=float, default=None)
    parser.add_argument(
        "--pressure-levels-hpa",
        default="1000,925,850,700,500",
        help="Comma-separated pressure levels for populated source-side atmospheres",
    )
    args = parser.parse_args()

    contract = SOURCE_CONTRACTS.get(args.source)
    if contract is None:
        raise SystemExit(
            f"Concrete prepared-case manifests are currently supported for "
            f"{', '.join(sorted(SOURCE_CONTRACTS))} only."
        )

    cycle_time = parse_utc(args.cycle_time_utc)
    offsets = build_offsets(contract.boundary_interval_seconds, args.forecast_hours)
    output_dir = (
        Path(args.output_dir)
        if args.output_dir is not None
        else REPO_ROOT / "cases" / "prepared" / f"{args.domain_name}_{args.source}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    tool_manifest_path = Path(args.manifest)
    if not tool_manifest_path.is_absolute():
        tool_manifest_path = REPO_ROOT / tool_manifest_path
    tool_manifest = parse_toolchain_manifest(tool_manifest_path)
    tool_records = make_tool_records(tool_manifest)

    grid = {
        "nx": args.nx,
        "ny": args.ny,
        "nz": args.nz,
        "dx": args.dx,
        "dy": args.dy,
        "z_top": args.z_top,
    }
    target_window = {
        "center_lat": args.center_lat,
        "center_lon": args.center_lon,
        "pressure_levels_hpa": parse_pressure_levels(args.pressure_levels_hpa),
    }
    times = {
        "cycle_time_utc": format_utc(cycle_time),
        "analysis_valid_time_utc": format_utc(cycle_time),
        "forecast_hours": args.forecast_hours,
    }

    analysis_stub = {
        "schema_version": "gwm-next-analysis-state/v1",
        "source": contract.display_name,
        "grid": grid,
        "cycle_time_utc": times["cycle_time_utc"],
        "valid_time_utc": times["analysis_valid_time_utc"],
        "forecast_offset_seconds": 0,
        "metadata": {
            "domain_name": args.domain_name,
            "prepared_at_utc": utc_now(),
            "status": "stub",
            "note": "Populate field values through external decode tooling before runtime ingest.",
            "target_window": target_window,
        },
        "field_groups": {
            "atmosphere": {
                "required_canonical_fields": canonical_names(contract.atmosphere),
                "field_bindings": binding_rows(contract.atmosphere),
                "storage": "external-population-required",
            },
            "surface": {
                "required_canonical_fields": canonical_names(contract.surface),
                "field_bindings": binding_rows(contract.surface),
                "storage": "external-population-required",
            },
            "static_surface": {
                "required_canonical_fields": canonical_names(contract.static_surface),
                "field_bindings": binding_rows(contract.static_surface),
                "storage": "external-population-required",
            },
        },
    }

    boundary_snapshots = []
    for offset in offsets:
        boundary_snapshots.append(
            {
                "forecast_offset_seconds": offset,
                "valid_time_utc": format_utc(cycle_time + timedelta(seconds=offset)),
                "field_groups": {
                    "atmosphere": {
                        "required_canonical_fields": canonical_names(contract.atmosphere),
                        "storage": "external-population-required",
                    },
                    "surface": {
                        "required_canonical_fields": canonical_names(contract.surface),
                        "storage": "external-population-required",
                    },
                },
            }
        )

    boundary_cache_stub = {
        "schema_version": "gwm-next-boundary-cache/v1",
        "source": contract.display_name,
        "grid": grid,
        "cycle_time_utc": times["cycle_time_utc"],
        "boundary_interval_seconds": contract.boundary_interval_seconds,
        "metadata": {
            "domain_name": args.domain_name,
            "prepared_at_utc": utc_now(),
            "status": "stub",
            "target_window": target_window,
        },
        "snapshots": boundary_snapshots,
    }

    checkpoint_stub = {
        "schema_version": "gwm-next-checkpoint-stub/v1",
        "restart_family": "gwm-native-checkpoint",
        "domain_name": args.domain_name,
        "source": contract.display_name,
        "cycle_time_utc": times["cycle_time_utc"],
        "cold_start_contract": {
            "analysis_state_schema": analysis_stub["schema_version"],
            "boundary_cache_schema": boundary_cache_stub["schema_version"],
            "required_groups": ["atmosphere", "surface", "static_surface"],
        },
        "runtime_products": {
            "verification_modes": ["screen", "aloft", "products"],
            "current_supported_bridge": "idealized_summary_json",
        },
        "note": "Checkpoint data is runtime-produced later; this file exists to lock the artifact contract now.",
    }

    product_plan = product_plan_stub(args.domain_name, contract)

    artifact_paths = {
        "analysis_state_stub": output_dir / "analysis_state_stub.json",
        "boundary_cache_stub": output_dir / "boundary_cache_stub.json",
        "checkpoint_stub": output_dir / "checkpoint_stub.json",
        "product_plan": output_dir / "product_plan.json",
        "prepared_case_manifest": output_dir / "prepared_case_manifest.json",
    }

    prepared_case_manifest = {
        "schema_version": "gwm-next-prepared-case/v1",
        "prepared_at_utc": utc_now(),
        "domain_name": args.domain_name,
        "source": {
            "key": contract.key,
            "display_name": contract.display_name,
            "first_class": contract.first_class,
            "boundary_interval_seconds": contract.boundary_interval_seconds,
        },
        "grid": grid,
        "target_window": target_window,
        "times": times,
        "contracts": {
            "analysis_state": {
                "schema_version": analysis_stub["schema_version"],
                "required_field_groups": {
                    "atmosphere": canonical_names(contract.atmosphere),
                    "surface": canonical_names(contract.surface),
                    "static_surface": canonical_names(contract.static_surface),
                },
            },
            "boundary_cache": {
                "schema_version": boundary_cache_stub["schema_version"],
                "interval_seconds": contract.boundary_interval_seconds,
                "forecast_offsets_seconds": offsets,
            },
            "checkpoint": {
                "schema_version": checkpoint_stub["schema_version"],
                "restart_family": checkpoint_stub["restart_family"],
            },
        },
        "external_toolchain": {
            "manifest_path": str(tool_manifest_path.resolve()),
            "manifest_exists": bool(tool_manifest.get("exists")),
            "policy": tool_manifest.get("policy", {}),
            "tools": tool_records,
        },
        "artifacts": {name: str(path.resolve()) for name, path in artifact_paths.items()},
    }

    for path, payload in (
        (artifact_paths["analysis_state_stub"], analysis_stub),
        (artifact_paths["boundary_cache_stub"], boundary_cache_stub),
        (artifact_paths["checkpoint_stub"], checkpoint_stub),
        (artifact_paths["product_plan"], product_plan),
        (artifact_paths["prepared_case_manifest"], prepared_case_manifest),
    ):
        path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")

    print(f"Prepared case manifest written to: {artifact_paths['prepared_case_manifest']}")
    print(f"Source contract: {contract.display_name}")
    print(f"Boundary snapshots: {len(boundary_snapshots)}")
    print("Artifacts:")
    for name, path in artifact_paths.items():
        print(f"- {name}: {path}")


if __name__ == "__main__":
    main()
