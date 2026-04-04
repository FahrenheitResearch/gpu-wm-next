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
        report = verify_idealized_summary(
            summary_path, load_json(summary_path)
        )
        child_reports.append(report)
    plan_view_relpath, plan_view_path = optional_members["plan_view"]
    if plan_view_path is not None:
        report = verify_plan_view_bundle(
            plan_view_path, load_json(plan_view_path)
        )
        child_reports.append(report)
    map_manifest_relpath, map_manifest_path = optional_members["map_manifest"]
    if map_manifest_path is not None:
        report = verify_map_manifest(
            map_manifest_path, load_json(map_manifest_path)
        )
        child_reports.append(report)

    rendered_field_consistency = True
    rendered_field_detail: dict[str, Any] = {"rendered": [], "plan_view_fields": []}
    if plan_view_path is not None and map_manifest_path is not None:
        plan_view_payload = load_json(plan_view_path)
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


if __name__ == "__main__":
    main()
