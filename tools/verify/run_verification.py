from __future__ import annotations

import argparse
import json
import math
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


BUNDLE_OPTIONAL_CANDIDATES = {
    "summary": ("summary.json", "runtime_summary.json"),
    "plan_view": ("plan_view.json", "runtime_plan_view.json"),
    "map_manifest": (
        "map_manifest.json",
        "plan_view_maps/map_manifest.json",
        "maps/map_manifest.json",
    ),
}

MOIST_FIELD_LIMITS = {
    "specific_humidity": {"min": 0.0, "max": 1.0},
    "specific_humidity_2m": {"min": 0.0, "max": 1.0},
    "cloud_water_mixing_ratio": {"min": 0.0, "max": 1.0},
    "rain_water_mixing_ratio": {"min": 0.0, "max": 1.0},
    "total_condensate": {"min": 0.0, "max": 1.0},
    "column_rain_water": {"min": 0.0, "max": 1.0e6},
    "accumulated_surface_precipitation": {"min": 0.0, "max": 1.0e6},
    "mean_surface_precipitation_rate": {"min": 0.0, "max": 1.0e6},
    "synthetic_reflectivity": {"min": -50.0, "max": 100.0},
    "relative_humidity": {"min": 0.0, "max": 150.0},
    "relative_humidity_2m": {"min": 0.0, "max": 150.0},
}


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def numeric_leaf_paths(value: Any, prefix: str = "") -> list[tuple[str, float]]:
    leaves: list[tuple[str, float]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            child_prefix = f"{prefix}.{key}" if prefix else str(key)
            leaves.extend(numeric_leaf_paths(child, child_prefix))
        return leaves
    if isinstance(value, list):
        for index, child in enumerate(value):
            child_prefix = f"{prefix}[{index}]"
            leaves.extend(numeric_leaf_paths(child, child_prefix))
        return leaves
    if is_finite_number(value):
        leaves.append((prefix or "<root>", float(value)))
    elif isinstance(value, (int, float)):
        leaves.append((prefix or "<root>", float("nan")))
    return leaves


def summary_tracer_path_sum(summary_payload: dict[str, Any], name: str) -> float | None:
    tracers = summary_payload.get("tracers", {})
    if not isinstance(tracers, dict):
        return None
    tracer_payload = tracers.get(name, {})
    if not isinstance(tracer_payload, dict):
        return None
    total_path = tracer_payload.get("column_integrated_sum_kg_m2")
    return float(total_path) if is_finite_number(total_path) else None


def summarize_plan_view_field(payload: dict[str, Any], field_name: str) -> dict[str, float] | None:
    for field in payload.get("fields", []):
        if str(field.get("name", "")) != field_name:
            continue
        values = field.get("values", [])
        if not all_finite(values):
            return None
        numeric_values = [float(value) for value in values]
        if not numeric_values:
            return {"sum": 0.0, "mean": 0.0, "max": 0.0}
        total = float(sum(numeric_values))
        return {
            "sum": total,
            "mean": total / float(len(numeric_values)),
            "max": max(numeric_values),
        }
    return None


def plan_view_field_values(payload: dict[str, Any], field_name: str) -> list[float] | None:
    for field in payload.get("fields", []):
        if str(field.get("name", "")) != field_name:
            continue
        values = field.get("values", [])
        if not all_finite(values):
            return None
        return [float(value) for value in values]
    return None


def all_finite(values: list[Any]) -> bool:
    return all(is_finite_number(value) for value in values)


def monotonic_offsets(offsets: list[int]) -> bool:
    return all(offsets[idx - 1] < offsets[idx] for idx in range(1, len(offsets)))


def list_bundle_files(path: Path) -> dict[str, Path]:
    return {
        str(child.relative_to(path)).replace("\\", "/"): child
        for child in path.rglob("*")
        if child.is_file()
    }


def find_bundle_member(files: dict[str, Path], candidates: tuple[str, ...]) -> tuple[str | None, Path | None]:
    for candidate in candidates:
        if candidate in files:
            return candidate, files[candidate]
    return None, None


def resolve_manifest_link(base_path: Path, raw_path: str | None) -> Path | None:
    if not raw_path:
        return None
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = (base_path.parent / candidate).resolve()
    return candidate


def verify_prepared_case_manifest(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    artifacts = payload.get("artifacts", {})
    contracts = payload.get("contracts", {})
    toolchain = payload.get("external_toolchain", {})

    required_group_counts = {
        "atmosphere": len(contracts.get("analysis_state", {}).get("required_field_groups", {}).get("atmosphere", [])),
        "surface": len(contracts.get("analysis_state", {}).get("required_field_groups", {}).get("surface", [])),
        "static_surface": len(contracts.get("analysis_state", {}).get("required_field_groups", {}).get("static_surface", [])),
    }
    offsets = contracts.get("boundary_cache", {}).get("forecast_offsets_seconds", [])
    checks = [
        {
            "name": "schema_version",
            "passed": payload.get("schema_version") == "gwm-next-prepared-case/v1",
            "detail": payload.get("schema_version", ""),
        },
        {
            "name": "required_field_groups_nonempty",
            "passed": all(count > 0 for count in required_group_counts.values()),
            "detail": required_group_counts,
        },
        {
            "name": "boundary_offsets_monotonic",
            "passed": isinstance(offsets, list) and bool(offsets) and monotonic_offsets(offsets),
            "detail": offsets,
        },
        {
            "name": "artifact_paths_exist",
            "passed": all(Path(value).exists() for value in artifacts.values()),
            "detail": artifacts,
        },
        {
            "name": "toolchain_manifest_present",
            "passed": bool(toolchain.get("manifest_exists")),
            "detail": toolchain.get("manifest_path", ""),
        },
    ]
    populated_artifact_keys = [
        name for name in ("analysis_state", "boundary_cache") if name in artifacts
    ]
    linked_reports: list[dict[str, Any]] = []
    if populated_artifact_keys:
        checks.append(
            {
                "name": "populated_artifact_pair_present",
                "passed": populated_artifact_keys == ["analysis_state", "boundary_cache"],
                "detail": populated_artifact_keys,
            }
        )
        for artifact_key, verifier in (
            ("analysis_state", verify_analysis_state_payload),
            ("boundary_cache", verify_boundary_cache_payload),
        ):
            artifact_path = artifacts.get(artifact_key)
            if artifact_path and Path(artifact_path).exists():
                linked_reports.append(
                    verifier(Path(artifact_path), load_json(Path(artifact_path)))
                )
        checks.append(
            {
                "name": "linked_runtime_artifacts_valid",
                "passed": all(report["overall_passed"] for report in linked_reports),
                "detail": [
                    {
                        "kind": report["input_kind"],
                        "overall_passed": report["overall_passed"],
                    }
                    for report in linked_reports
                ],
            }
        )
    tool_records = toolchain.get("tools", [])
    return {
        "input_kind": "prepared_case_manifest",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "linked_reports": linked_reports,
        "available_external_tools": [
            record["name"] for record in tool_records if record.get("exists")
        ],
        "missing_external_tools": [
            record["name"] for record in tool_records if not record.get("exists")
        ],
        "overall_passed": all(check["passed"] for check in checks),
    }


def _required_field_group_checks(
    values: dict[str, list[Any]],
    expected_size: int,
    groups: list[tuple[str, ...]],
    label: str,
) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    for group in groups:
        present = [name for name in group if name in values]
        sizes = {name: len(values[name]) for name in present}
        passed = bool(present) and all(size == expected_size for size in sizes.values())
        checks.append(
            {
                "name": f"{label}_{'_or_'.join(group)}",
                "passed": passed,
                "detail": {
                    "required_any_of": list(group),
                    "present": present,
                    "expected_size": expected_size,
                    "sizes": sizes,
                },
            }
        )
    return checks


def verify_analysis_state_payload(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    grid = payload.get("grid", {})
    metadata = payload.get("metadata", {})
    atmosphere = payload.get("atmosphere", {})
    surface = payload.get("surface", {})
    static_surface = payload.get("static_surface", {})
    expected_3d = int(grid.get("nx", 0)) * int(grid.get("ny", 0)) * int(grid.get("nz", 0))
    expected_2d = int(grid.get("nx", 0)) * int(grid.get("ny", 0))
    checks = [
        {
            "name": "schema_version",
            "passed": payload.get("schema_version") == "gwm-next-analysis-state/v1",
            "detail": payload.get("schema_version", ""),
        },
        {
            "name": "population_status",
            "passed": metadata.get("status") == "populated",
            "detail": metadata.get("status", ""),
        },
        {
            "name": "grid_shape_positive",
            "passed": int(grid.get("nx", 0)) > 0
            and int(grid.get("ny", 0)) > 0
            and int(grid.get("nz", 0)) > 0
            and float(grid.get("z_top", 0.0)) > 0.0,
            "detail": grid,
        },
        {
            "name": "source_present",
            "passed": bool(payload.get("source")),
            "detail": payload.get("source", ""),
        },
        {
            "name": "target_window_present",
            "passed": bool(metadata.get("target_window")),
            "detail": metadata.get("target_window", {}),
        },
    ]
    checks.extend(
        _required_field_group_checks(
            atmosphere,
            expected_3d,
            [
                ("u_wind",),
                ("v_wind",),
                ("w_wind",),
                ("air_temperature",),
                ("air_pressure",),
                ("geopotential_height",),
                ("specific_humidity", "water_vapor_mixing_ratio"),
            ],
            "atmosphere",
        )
    )
    checks.extend(
        _required_field_group_checks(
            surface,
            expected_2d,
            [
                ("surface_pressure",),
                ("air_temperature_2m",),
                ("specific_humidity_2m",),
                ("u_wind_10m",),
                ("v_wind_10m",),
                ("skin_temperature",),
            ],
            "surface",
        )
    )
    checks.extend(
        _required_field_group_checks(
            static_surface,
            expected_2d,
            [("terrain_height",), ("land_mask",), ("land_use_index",)],
            "static_surface",
        )
    )
    return {
        "input_kind": "analysis_state",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "field_counts": {
            "atmosphere": len(atmosphere),
            "surface": len(surface),
            "static_surface": len(static_surface),
        },
        "overall_passed": all(check["passed"] for check in checks),
    }


def verify_boundary_cache_payload(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    grid = payload.get("grid", {})
    metadata = payload.get("metadata", {})
    snapshots = payload.get("snapshots", [])
    expected_3d = int(grid.get("nx", 0)) * int(grid.get("ny", 0)) * int(grid.get("nz", 0))
    expected_2d = int(grid.get("nx", 0)) * int(grid.get("ny", 0))
    offsets = [snapshot.get("forecast_offset_seconds") for snapshot in snapshots]
    checks = [
        {
            "name": "schema_version",
            "passed": payload.get("schema_version") == "gwm-next-boundary-cache/v1",
            "detail": payload.get("schema_version", ""),
        },
        {
            "name": "population_status",
            "passed": metadata.get("status") == "populated",
            "detail": metadata.get("status", ""),
        },
        {
            "name": "grid_shape_positive",
            "passed": int(grid.get("nx", 0)) > 0
            and int(grid.get("ny", 0)) > 0
            and int(grid.get("nz", 0)) > 0,
            "detail": grid,
        },
        {
            "name": "boundary_interval_positive",
            "passed": int(payload.get("boundary_interval_seconds", 0)) > 0,
            "detail": payload.get("boundary_interval_seconds", 0),
        },
        {
            "name": "source_present",
            "passed": bool(payload.get("source")),
            "detail": payload.get("source", ""),
        },
        {
            "name": "target_window_present",
            "passed": bool(metadata.get("target_window")),
            "detail": metadata.get("target_window", {}),
        },
        {
            "name": "snapshots_present",
            "passed": len(snapshots) >= 2,
            "detail": len(snapshots),
        },
        {
            "name": "offsets_monotonic",
            "passed": all(isinstance(offset, int) for offset in offsets) and monotonic_offsets(offsets),
            "detail": offsets,
        },
    ]
    for snapshot_index, snapshot in enumerate(snapshots):
        atmosphere = snapshot.get("atmosphere", {})
        surface = snapshot.get("surface", {})
        checks.extend(
            _required_field_group_checks(
                atmosphere,
                expected_3d,
                [
                    ("u_wind",),
                    ("v_wind",),
                    ("w_wind",),
                    ("air_temperature",),
                    ("air_pressure",),
                    ("geopotential_height",),
                    ("specific_humidity", "water_vapor_mixing_ratio"),
                ],
                f"snapshot_{snapshot_index}_atmosphere",
            )
        )
        checks.extend(
            _required_field_group_checks(
                surface,
                expected_2d,
                [
                    ("surface_pressure",),
                    ("air_temperature_2m",),
                    ("specific_humidity_2m",),
                    ("u_wind_10m",),
                    ("v_wind_10m",),
                    ("skin_temperature",),
                ],
                f"snapshot_{snapshot_index}_surface",
            )
        )
        checks.append(
            {
                "name": f"snapshot_{snapshot_index}_valid_time_present",
                "passed": bool(snapshot.get("valid_time_utc")),
                "detail": snapshot.get("valid_time_utc", ""),
            }
        )
    return {
        "input_kind": "boundary_cache",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "snapshot_count": len(snapshots),
        "overall_passed": all(check["passed"] for check in checks),
    }


def verify_source_run_bundle(path: Path) -> dict[str, Any]:
    members = list_bundle_files(path)
    child_reports: list[dict[str, Any]] = []
    summary_payload: dict[str, Any] | None = None
    summary_kind: str | None = None
    plan_view_payload: dict[str, Any] | None = None

    required_members = [
        "prepared_case_manifest.json",
        "analysis_state.json",
        "boundary_cache.json",
    ]
    optional_members = {
        name: find_bundle_member(members, candidates)
        for name, candidates in BUNDLE_OPTIONAL_CANDIDATES.items()
    }

    checks = [
        {
            "name": "bundle_is_directory",
            "passed": path.is_dir(),
            "detail": str(path.resolve()),
        },
        {
            "name": "required_members_present",
            "passed": all(member in members for member in required_members),
            "detail": {member: (member in members) for member in required_members},
        },
    ]

    if "prepared_case_manifest.json" in members:
        report = verify_prepared_case_manifest(
            members["prepared_case_manifest.json"],
            load_json(members["prepared_case_manifest.json"]),
        )
        child_reports.append(report)
    if "analysis_state.json" in members:
        report = verify_analysis_state_payload(
            members["analysis_state.json"], load_json(members["analysis_state.json"])
        )
        child_reports.append(report)
    if "boundary_cache.json" in members:
        report = verify_boundary_cache_payload(
            members["boundary_cache.json"], load_json(members["boundary_cache.json"])
        )
        child_reports.append(report)
    summary_relpath, summary_path = optional_members["summary"]
    if summary_path is not None:
        summary_payload = load_json(summary_path)
        summary_kind = detect_kind(summary_payload)
        if summary_kind == "runtime_summary":
            report = verify_runtime_summary(summary_path, summary_payload)
        else:
            report = verify_idealized_summary(summary_path, summary_payload)
        child_reports.append(report)
    plan_view_relpath, plan_view_path = optional_members["plan_view"]
    if plan_view_path is not None:
        plan_view_payload = load_json(plan_view_path)
        report = verify_plan_view_bundle(plan_view_path, plan_view_payload)
        child_reports.append(report)
    map_manifest_relpath, map_manifest_path = optional_members["map_manifest"]
    if map_manifest_path is not None:
        report = verify_map_manifest(
            map_manifest_path, load_json(map_manifest_path)
        )
        child_reports.append(report)

    rendered_field_consistency = True
    rendered_field_detail: dict[str, Any] = {"rendered": [], "plan_view_fields": []}
    if plan_view_path is not None and map_manifest_path is not None and plan_view_payload is not None:
        map_manifest_payload = load_json(map_manifest_path)
        plan_view_fields = {
            str(field.get("name", ""))
            for field in plan_view_payload.get("fields", [])
            if field.get("name")
        }
        rendered_fields = {
            str(image.get("source_field_name") or image.get("field", ""))
            for image in map_manifest_payload.get("images", [])
            if image.get("source_field_name") or image.get("field")
        }
        rendered_field_consistency = rendered_fields.issubset(plan_view_fields)
        rendered_field_detail = {
            "rendered": sorted(rendered_fields),
            "plan_view_fields": sorted(plan_view_fields),
            "plan_view_path": plan_view_relpath,
            "map_manifest_path": map_manifest_relpath,
        }

    summary_plan_view_consistent = True
    summary_plan_view_detail: dict[str, Any] = {}
    if summary_payload is not None and plan_view_payload is not None:
        if summary_kind == "runtime_summary":
            final_moisture = summary_payload.get("final", {}).get("moisture", {})
            accumulated_field = summarize_plan_view_field(
                plan_view_payload, "accumulated_surface_precipitation"
            )
            rate_field = summarize_plan_view_field(
                plan_view_payload, "mean_surface_precipitation_rate"
            )
            accumulated_values = plan_view_field_values(
                plan_view_payload, "accumulated_surface_precipitation"
            )
            rate_values = plan_view_field_values(
                plan_view_payload, "mean_surface_precipitation_rate"
            )
            plan_view_fields = {
                str(field.get("name", ""))
                for field in plan_view_payload.get("fields", [])
                if field.get("name")
            }
            elapsed_hours = float(summary_payload.get("elapsed_hours", 0.0))
            if elapsed_hours <= 0.0:
                elapsed_hours = (
                    float(summary_payload.get("steps", 0))
                    * float(summary_payload.get("dt", 0.0))
                ) / 3600.0
            summary_plan_view_consistent = (
                isinstance(final_moisture, dict) and accumulated_field is not None
            )
            if summary_plan_view_consistent and isinstance(final_moisture, dict):
                summary_plan_view_consistent = (
                    abs(
                        accumulated_field["sum"]
                        - float(
                            final_moisture.get(
                                "accumulated_surface_precipitation_sum_mm", 0.0
                            )
                        )
                    )
                    <= 1.0e-6
                    and abs(
                        accumulated_field["mean"]
                        - float(
                            final_moisture.get(
                                "mean_surface_precipitation_mm", 0.0
                            )
                        )
                    )
                    <= 1.0e-6
                    and abs(
                        accumulated_field["max"]
                        - float(
                            final_moisture.get(
                                "max_surface_precipitation_mm", 0.0
                            )
                        )
                    )
                    <= 1.0e-6
                )
                if rate_field is not None and elapsed_hours > 0.0:
                    summary_plan_view_consistent = summary_plan_view_consistent and (
                        abs(
                            rate_field["mean"]
                            - float(
                                final_moisture.get(
                                    "mean_surface_precipitation_mm", 0.0
                                )
                            )
                            / elapsed_hours
                        )
                        <= 1.0e-6
                    )
                    if accumulated_values is not None and rate_values is not None:
                        summary_plan_view_consistent = (
                            summary_plan_view_consistent
                            and len(accumulated_values) == len(rate_values)
                            and all(
                                abs(rate - accum / elapsed_hours) <= 1.0e-6
                                for accum, rate in zip(accumulated_values, rate_values)
                            )
                        )
                if float(final_moisture.get("condensed_water_path_sum_kg_m2", 0.0)) > 0.0:
                    summary_plan_view_consistent = summary_plan_view_consistent and all(
                        field_name in plan_view_fields
                        for field_name in (
                            "cloud_water_mixing_ratio",
                            "rain_water_mixing_ratio",
                            "total_condensate",
                            "column_cloud_water",
                            "column_rain_water",
                            "column_total_condensate",
                            "column_rain_fraction",
                            "synthetic_reflectivity",
                            "accumulated_surface_precipitation",
                            "mean_surface_precipitation_rate",
                        )
                    )
            summary_plan_view_detail = {
                "elapsed_hours": elapsed_hours,
                "plan_view_fields": sorted(plan_view_fields),
                "accumulated_surface_precipitation": accumulated_field,
                "mean_surface_precipitation_rate": rate_field,
                "cellwise_rate_matches_accumulation": (
                    accumulated_values is not None
                    and rate_values is not None
                    and elapsed_hours > 0.0
                    and len(accumulated_values) == len(rate_values)
                    and all(
                        abs(rate - accum / elapsed_hours) <= 1.0e-6
                        for accum, rate in zip(accumulated_values, rate_values)
                    )
                ),
                "final_moisture": final_moisture,
            }

    checks.extend(
        [
            {
                "name": "optional_members_tracked",
                "passed": True,
                "detail": {
                    member: relpath
                    for member, (relpath, _) in optional_members.items()
                    if relpath is not None
                },
            },
            {
                "name": "child_reports_passed",
                "passed": all(report["overall_passed"] for report in child_reports),
                "detail": [
                    {
                        "kind": report["input_kind"],
                        "overall_passed": report["overall_passed"],
                    }
                    for report in child_reports
                ],
            },
            {
                "name": "rendered_fields_match_plan_view",
                "passed": rendered_field_consistency,
                "detail": rendered_field_detail,
            },
            {
                "name": "runtime_summary_matches_plan_view",
                "passed": summary_plan_view_consistent,
                "detail": summary_plan_view_detail,
            },
        ]
    )

    return {
        "input_kind": "source_run_bundle",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "bundle_members": sorted(members),
        "child_reports": child_reports,
        "overall_passed": all(check["passed"] for check in checks),
    }


def verify_boundary_cache_stub(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    snapshots = payload.get("snapshots", [])
    offsets = [snapshot.get("forecast_offset_seconds") for snapshot in snapshots]
    checks = [
        {
            "name": "schema_version",
            "passed": payload.get("schema_version") == "gwm-next-boundary-cache/v1",
            "detail": payload.get("schema_version", ""),
        },
        {
            "name": "snapshots_present",
            "passed": len(snapshots) >= 2,
            "detail": len(snapshots),
        },
        {
            "name": "offsets_monotonic",
            "passed": all(isinstance(offset, int) for offset in offsets) and monotonic_offsets(offsets),
            "detail": offsets,
        },
    ]
    return {
        "input_kind": "boundary_cache_stub",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "overall_passed": all(check["passed"] for check in checks),
    }


def verify_checkpoint_stub(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    cold_start = payload.get("cold_start_contract", {})
    checks = [
        {
            "name": "schema_version",
            "passed": payload.get("schema_version") == "gwm-next-checkpoint-stub/v1",
            "detail": payload.get("schema_version", ""),
        },
        {
            "name": "restart_family",
            "passed": payload.get("restart_family") == "gwm-native-checkpoint",
            "detail": payload.get("restart_family", ""),
        },
        {
            "name": "cold_start_groups",
            "passed": cold_start.get("required_groups") == ["atmosphere", "surface", "static_surface"],
            "detail": cold_start.get("required_groups", []),
        },
    ]
    return {
        "input_kind": "checkpoint_stub",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "overall_passed": all(check["passed"] for check in checks),
    }


def verify_product_plan(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    modes = payload.get("verification_modes", [])
    planned = payload.get("planned_products", {})
    checks = [
        {
            "name": "schema_version",
            "passed": payload.get("schema_version") == "gwm-next-product-plan/v1",
            "detail": payload.get("schema_version", ""),
        },
        {
            "name": "verification_modes",
            "passed": modes == ["screen", "aloft", "products"],
            "detail": modes,
        },
        {
            "name": "planned_products_nonempty",
            "passed": all(bool(planned.get(section)) for section in ("screen", "aloft", "products")),
            "detail": planned,
        },
    ]
    return {
        "input_kind": "product_plan",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "overall_passed": all(check["passed"] for check in checks),
    }


def verify_map_manifest(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    images = payload.get("images", [])
    image_fields = [
        str(image.get("source_field_name") or image.get("field", ""))
        for image in images
        if image.get("source_field_name") or image.get("field")
    ]
    plan_view_path = resolve_manifest_link(path, payload.get("input_plan_view_path"))
    plan_view_exists = bool(plan_view_path and plan_view_path.exists())
    plan_view_fields: list[str] = []
    plan_view_subset_passed = True
    if plan_view_exists and plan_view_path is not None:
        plan_view_payload = load_json(plan_view_path)
        plan_view_fields = [
            str(field.get("name", ""))
            for field in plan_view_payload.get("fields", [])
            if field.get("name")
        ]
        plan_view_subset_passed = set(image_fields).issubset(set(plan_view_fields))
    checks = [
        {
            "name": "schema_version",
            "passed": payload.get("schema_version") == "gwm-next-map-manifest/v1",
            "detail": payload.get("schema_version", ""),
        },
        {
            "name": "images_present",
            "passed": bool(images),
            "detail": len(images),
        },
        {
            "name": "image_paths_exist",
            "passed": all(Path(image.get("path", "")).exists() for image in images),
            "detail": [image.get("path", "") for image in images],
        },
        {
            "name": "image_fields_unique",
            "passed": len(image_fields) == len(set(image_fields)),
            "detail": image_fields,
        },
        {
            "name": "input_plan_view_exists",
            "passed": plan_view_exists,
            "detail": str(plan_view_path.resolve()) if plan_view_exists and plan_view_path else str(payload.get("input_plan_view_path", "")),
        },
        {
            "name": "rendered_fields_subset_of_plan_view",
            "passed": plan_view_subset_passed,
            "detail": {
                "rendered_fields": image_fields,
                "plan_view_fields": plan_view_fields,
            },
        },
    ]
    return {
        "input_kind": "map_manifest",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "overall_passed": all(check["passed"] for check in checks),
    }


def verify_plan_view_bundle(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    fields = payload.get("fields", [])
    grid = payload.get("grid", {})
    field_names = [str(field.get("name", "")) for field in fields if field.get("name")]
    finite_value_checks = []
    moist_field_details: list[dict[str, Any]] = []
    moist_field_bounds_passed = True
    for field in fields:
        values = field.get("values", [])
        finite_value_checks.append(all_finite(values))
        name = str(field.get("name", ""))
        if name in MOIST_FIELD_LIMITS:
            numeric_values = [float(value) for value in values]
            if numeric_values:
                lower = min(numeric_values)
                upper = max(numeric_values)
            else:
                lower = 0.0
                upper = 0.0
            limits = MOIST_FIELD_LIMITS[name]
            field_passed = lower >= limits["min"] - 1.0e-6 and upper <= limits["max"] + 1.0e-6
            moist_field_bounds_passed = moist_field_bounds_passed and field_passed
            moist_field_details.append(
                {
                    "name": name,
                    "min": lower,
                    "max": upper,
                    "units": field.get("units", ""),
                    "passed": field_passed,
                }
            )
    checks = [
        {
            "name": "schema_version",
            "passed": payload.get("schema_version") == "gwm-next-plan-view/v1",
            "detail": payload.get("schema_version", ""),
        },
        {
            "name": "fields_present",
            "passed": bool(fields),
            "detail": len(fields),
        },
        {
            "name": "grid_shape_present",
            "passed": int(grid.get("nx", 0)) > 0 and int(grid.get("ny", 0)) > 0,
            "detail": grid,
        },
        {
            "name": "field_storage_consistent",
            "passed": all(
                int(field.get("nx", 0)) * int(field.get("ny", 0)) == len(field.get("values", []))
                for field in fields
            ),
            "detail": [field.get("name", "") for field in fields],
        },
        {
            "name": "field_names_unique",
            "passed": len(field_names) == len(set(field_names)),
            "detail": field_names,
        },
        {
            "name": "field_values_finite",
            "passed": all(finite_value_checks),
            "detail": field_names,
        },
        {
            "name": "moist_field_bounds",
            "passed": moist_field_bounds_passed,
            "detail": moist_field_details,
        },
    ]
    return {
        "input_kind": "plan_view_bundle",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "field_names": field_names,
        "moist_fields_present": [detail["name"] for detail in moist_field_details],
        "overall_passed": all(check["passed"] for check in checks),
    }


def verify_runtime_summary(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    initial = payload.get("initial", {})
    final = payload.get("final", {})
    initial_numeric_leaves = numeric_leaf_paths(initial, "initial")
    final_numeric_leaves = numeric_leaf_paths(final, "final")
    elapsed_seconds = payload.get("elapsed_seconds")
    if not is_finite_number(elapsed_seconds):
        elapsed_seconds = float(payload.get("steps", 0)) * float(payload.get("dt", 0.0))

    moisture_fields = (
        "vapor_water_path_sum_kg_m2",
        "cloud_water_path_sum_kg_m2",
        "rain_water_path_sum_kg_m2",
        "condensed_water_path_sum_kg_m2",
        "total_water_path_sum_kg_m2",
        "accumulated_surface_precipitation_sum_mm",
        "mean_surface_precipitation_mm",
        "max_surface_precipitation_mm",
        "total_surface_cell_count",
        "precipitating_surface_cell_count",
        "precipitating_surface_fraction",
        "mean_precipitating_surface_precipitation_mm",
    )
    initial_moisture = initial.get("moisture", {})
    final_moisture = final.get("moisture", {})
    tracer_names = (
        "specific_humidity",
        "cloud_water_mixing_ratio",
        "rain_water_mixing_ratio",
    )

    moisture_details: list[dict[str, Any]] = []
    moisture_fields_present = True
    tracer_consistency = True
    for label, summary_payload in (("initial", initial), ("final", final)):
        moisture = summary_payload.get("moisture", {})
        moisture_fields_present = moisture_fields_present and isinstance(moisture, dict)
        if isinstance(moisture, dict):
            moisture_details.append(
                {
                    "label": label,
                    "missing_fields": [
                        field for field in moisture_fields if field not in moisture
                    ],
                }
            )
        else:
            moisture_details.append({"label": label, "missing_fields": list(moisture_fields)})
        tracer_qv = summary_tracer_path_sum(summary_payload, "specific_humidity")
        tracer_qc = summary_tracer_path_sum(summary_payload, "cloud_water_mixing_ratio")
        tracer_qr = summary_tracer_path_sum(summary_payload, "rain_water_mixing_ratio")
        if isinstance(moisture, dict):
            if tracer_qv is not None and is_finite_number(moisture.get("vapor_water_path_sum_kg_m2")):
                tracer_consistency = tracer_consistency and (
                    abs(float(moisture["vapor_water_path_sum_kg_m2"]) - tracer_qv) <= 1.0e-6
                )
            if tracer_qc is not None and is_finite_number(moisture.get("cloud_water_path_sum_kg_m2")):
                tracer_consistency = tracer_consistency and (
                    abs(float(moisture["cloud_water_path_sum_kg_m2"]) - tracer_qc) <= 1.0e-6
                )
            if tracer_qr is not None and is_finite_number(moisture.get("rain_water_path_sum_kg_m2")):
                tracer_consistency = tracer_consistency and (
                    abs(float(moisture["rain_water_path_sum_kg_m2"]) - tracer_qr) <= 1.0e-6
                )

    water_path_budget_consistent = False
    water_path_detail: dict[str, Any] = {}
    if isinstance(final_moisture, dict) and all(
        key in final_moisture for key in moisture_fields[:5]
    ):
        vapor_path = float(final_moisture["vapor_water_path_sum_kg_m2"])
        cloud_path = float(final_moisture["cloud_water_path_sum_kg_m2"])
        rain_path = float(final_moisture["rain_water_path_sum_kg_m2"])
        condensed_path = float(final_moisture["condensed_water_path_sum_kg_m2"])
        total_path = float(final_moisture["total_water_path_sum_kg_m2"])
        water_path_budget_consistent = (
            abs(condensed_path - (cloud_path + rain_path)) <= 1.0e-6
            and abs(total_path - (vapor_path + condensed_path)) <= 1.0e-6
        )
        water_path_detail = {
            "vapor_path": vapor_path,
            "cloud_path": cloud_path,
            "rain_path": rain_path,
            "condensed_path": condensed_path,
            "total_path": total_path,
        }

    precip_summary_consistent = False
    precip_detail: dict[str, Any] = {}
    if isinstance(final_moisture, dict) and all(
        key in final_moisture for key in moisture_fields
    ):
        total_cells = int(final_moisture["total_surface_cell_count"])
        precip_cells = int(final_moisture["precipitating_surface_cell_count"])
        precip_sum = float(final_moisture["accumulated_surface_precipitation_sum_mm"])
        mean_precip = float(final_moisture["mean_surface_precipitation_mm"])
        max_precip = float(final_moisture["max_surface_precipitation_mm"])
        precip_fraction = float(final_moisture["precipitating_surface_fraction"])
        wet_mean_precip = float(
            final_moisture["mean_precipitating_surface_precipitation_mm"]
        )
        precip_summary_consistent = (
            total_cells >= 0
            and 0 <= precip_cells <= total_cells
            and precip_sum >= -1.0e-9
            and mean_precip >= -1.0e-9
            and max_precip >= -1.0e-9
            and wet_mean_precip >= -1.0e-9
            and abs(mean_precip * float(total_cells) - precip_sum) <= 1.0e-6
            and (
                total_cells == 0
                or abs(precip_fraction - float(precip_cells) / float(total_cells))
                <= 1.0e-6
            )
            and (
                precip_cells == 0
                or abs(wet_mean_precip * float(precip_cells) - precip_sum) <= 1.0e-6
            )
            and (precip_cells > 0 or abs(wet_mean_precip) <= 1.0e-9)
            and max_precip + 1.0e-9 >= wet_mean_precip
            and wet_mean_precip + 1.0e-9 >= mean_precip
        )
        precip_detail = {
            "elapsed_seconds": elapsed_seconds,
            "precip_sum_mm": precip_sum,
            "mean_precip_mm": mean_precip,
            "max_precip_mm": max_precip,
            "total_cells": total_cells,
            "precipitating_cells": precip_cells,
            "precipitating_fraction": precip_fraction,
            "wet_mean_precip_mm": wet_mean_precip,
        }

    accumulated_monotonic = False
    if isinstance(initial_moisture, dict) and isinstance(final_moisture, dict):
        if is_finite_number(initial_moisture.get("accumulated_surface_precipitation_sum_mm")) and is_finite_number(
            final_moisture.get("accumulated_surface_precipitation_sum_mm")
        ):
            accumulated_monotonic = (
                float(final_moisture["accumulated_surface_precipitation_sum_mm"])
                + 1.0e-9
                >= float(initial_moisture["accumulated_surface_precipitation_sum_mm"])
            )

    checks = [
        {
            "name": "runtime_summary_shape",
            "passed": all(
                key in payload
                for key in (
                    "case",
                    "source",
                    "cycle_time_utc",
                    "analysis_valid_time_utc",
                    "steps",
                    "dt",
                    "fast_substeps",
                    "grid",
                    "surface_runtime",
                    "initial",
                    "final",
                )
            ),
            "detail": sorted(payload.keys()),
        },
        {
            "name": "initial_numeric_leaves_finite",
            "passed": all(math.isfinite(value) for _, value in initial_numeric_leaves),
            "detail": [leaf_path for leaf_path, _ in initial_numeric_leaves],
        },
        {
            "name": "final_numeric_leaves_finite",
            "passed": all(math.isfinite(value) for _, value in final_numeric_leaves),
            "detail": [leaf_path for leaf_path, _ in final_numeric_leaves],
        },
        {
            "name": "moisture_fields_present",
            "passed": moisture_fields_present
            and all(not detail["missing_fields"] for detail in moisture_details),
            "detail": moisture_details,
        },
        {
            "name": "tracer_moisture_consistent",
            "passed": tracer_consistency,
            "detail": {"checked_tracers": list(tracer_names)},
        },
        {
            "name": "water_path_budget_consistent",
            "passed": water_path_budget_consistent,
            "detail": water_path_detail,
        },
        {
            "name": "precipitation_summary_consistent",
            "passed": precip_summary_consistent,
            "detail": precip_detail,
        },
        {
            "name": "accumulated_precipitation_monotonic",
            "passed": accumulated_monotonic,
            "detail": {
                "initial": initial_moisture.get("accumulated_surface_precipitation_sum_mm"),
                "final": final_moisture.get("accumulated_surface_precipitation_sum_mm"),
            },
        },
    ]

    return {
        "input_kind": "runtime_summary",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "overall_passed": all(check["passed"] for check in checks),
    }


def verify_idealized_summary(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    initial = payload.get("initial", {})
    final = payload.get("final", {})
    initial_numeric_leaves = numeric_leaf_paths(initial, "initial")
    final_numeric_leaves = numeric_leaf_paths(final, "final")
    total_rho_initial = initial.get("total_dry_mass")
    total_rho_final = final.get("total_dry_mass")
    total_theta_initial = initial.get("total_rho_theta_m")
    total_theta_final = final.get("total_rho_theta_m")
    w_face = final.get("w_face", {})
    checks = [
        {
            "name": "summary_shape",
            "passed": all(key in payload for key in ("case", "steps", "dt", "initial", "final")),
            "detail": sorted(payload.keys()),
        },
        {
            "name": "mass_totals_finite",
            "passed": is_finite_number(total_rho_initial) and is_finite_number(total_rho_final),
            "detail": {"initial": total_rho_initial, "final": total_rho_final},
        },
        {
            "name": "theta_totals_finite",
            "passed": is_finite_number(total_theta_initial) and is_finite_number(total_theta_final),
            "detail": {"initial": total_theta_initial, "final": total_theta_final},
        },
        {
            "name": "w_face_finite",
            "passed": is_finite_number(w_face.get("min")) and is_finite_number(w_face.get("max")),
            "detail": w_face,
        },
        {
            "name": "initial_numeric_leaves_finite",
            "passed": all(math.isfinite(value) for _, value in initial_numeric_leaves),
            "detail": [path for path, _ in initial_numeric_leaves],
        },
        {
            "name": "final_numeric_leaves_finite",
            "passed": all(math.isfinite(value) for _, value in final_numeric_leaves),
            "detail": [path for path, _ in final_numeric_leaves],
        },
    ]
    return {
        "input_kind": "idealized_summary",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "derived_metrics": {
            "mass_drift": None
            if not all(is_finite_number(v) for v in (total_rho_initial, total_rho_final))
            else float(total_rho_final) - float(total_rho_initial),
            "rho_theta_drift": None
            if not all(is_finite_number(v) for v in (total_theta_initial, total_theta_final))
            else float(total_theta_final) - float(total_theta_initial),
            "w_face_span": None
            if not (is_finite_number(w_face.get("min")) and is_finite_number(w_face.get("max")))
            else float(w_face["max"]) - float(w_face["min"]),
        },
        "product_bridge": {
            "current_supported_input": "idealized_summary_json",
            "verification_modes": ["products"],
            "available_products": ["dry_state_summary_report"],
            "planned_follow_on_products": ["map_manifest_json", "screen_panel_json"],
        },
        "overall_passed": all(check["passed"] for check in checks),
    }


def detect_kind(payload: dict[str, Any]) -> str:
    schema_version = payload.get("schema_version")
    if schema_version == "gwm-next-prepared-case/v1":
        return "prepared_case_manifest"
    if schema_version == "gwm-next-analysis-state/v1":
        return "analysis_state"
    if schema_version == "gwm-next-boundary-cache/v1":
        metadata = payload.get("metadata", {})
        if metadata.get("status") == "stub":
            return "boundary_cache_stub"
        return "boundary_cache"
    if schema_version == "gwm-next-checkpoint-stub/v1":
        return "checkpoint_stub"
    if schema_version == "gwm-next-product-plan/v1":
        return "product_plan"
    if schema_version == "gwm-next-plan-view/v1":
        return "plan_view_bundle"
    if schema_version == "gwm-next-map-manifest/v1":
        return "map_manifest"
    if all(key in payload for key in ("case", "initial", "final")):
        if all(
            key in payload
            for key in (
                "source",
                "cycle_time_utc",
                "analysis_valid_time_utc",
                "grid",
                "surface_runtime",
            )
        ):
            return "runtime_summary"
        return "idealized_summary"
    raise ValueError("Unable to detect input kind from JSON payload")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Verify prepared-case manifests, populated analysis/boundary artifacts, "
            "source-run bundles, checkpoint stubs, and idealized summary products."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="JSON artifact path or prepared/source-run bundle directory to verify",
    )
    parser.add_argument(
        "--kind",
        default="auto",
        choices=[
            "auto",
            "prepared_case_manifest",
            "analysis_state",
            "boundary_cache",
            "boundary_cache_stub",
            "checkpoint_stub",
            "product_plan",
            "plan_view_bundle",
            "map_manifest",
            "runtime_summary",
            "idealized_summary",
            "source_run_bundle",
        ],
    )
    parser.add_argument("--output-json", default=None, help="Optional explicit report path")
    args = parser.parse_args()

    input_path = Path(args.input)
    if input_path.is_dir():
        kind = "source_run_bundle" if args.kind == "auto" else args.kind
        if kind != "source_run_bundle":
            raise ValueError("Directory inputs only support kind=auto or source_run_bundle")
        report = verify_source_run_bundle(input_path)
    else:
        payload = load_json(input_path)
        kind = detect_kind(payload) if args.kind == "auto" else args.kind

        if kind == "prepared_case_manifest":
            report = verify_prepared_case_manifest(input_path, payload)
        elif kind == "analysis_state":
            report = verify_analysis_state_payload(input_path, payload)
        elif kind == "boundary_cache":
            report = verify_boundary_cache_payload(input_path, payload)
        elif kind == "boundary_cache_stub":
            report = verify_boundary_cache_stub(input_path, payload)
        elif kind == "checkpoint_stub":
            report = verify_checkpoint_stub(input_path, payload)
        elif kind == "product_plan":
            report = verify_product_plan(input_path, payload)
        elif kind == "plan_view_bundle":
            report = verify_plan_view_bundle(input_path, payload)
        elif kind == "map_manifest":
            report = verify_map_manifest(input_path, payload)
        elif kind == "runtime_summary":
            report = verify_runtime_summary(input_path, payload)
        elif kind == "idealized_summary":
            report = verify_idealized_summary(input_path, payload)
        else:  # pragma: no cover
            raise ValueError(f"Unsupported verification kind: {kind}")

    output_path = (
        Path(args.output_json)
        if args.output_json is not None
        else (input_path / "verification.json" if input_path.is_dir() else input_path.with_suffix(".verification.json"))
    )
    output_path.write_text(json.dumps(report, indent=2, sort_keys=False) + "\n", encoding="utf-8")

    print(f"Verification kind: {kind}")
    print(f"Input: {input_path.resolve()}")
    print(f"Report: {output_path.resolve()}")
    print(f"Overall passed: {report['overall_passed']}")
    for check in report["checks"]:
        status = "PASS" if check["passed"] else "FAIL"
        print(f"- [{status}] {check['name']}")
    if not report["overall_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
