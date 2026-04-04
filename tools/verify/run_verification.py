from __future__ import annotations

import argparse
import json
import math
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def monotonic_offsets(offsets: list[int]) -> bool:
    return all(offsets[idx - 1] < offsets[idx] for idx in range(1, len(offsets)))


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
    tool_records = toolchain.get("tools", [])
    return {
        "input_kind": "prepared_case_manifest",
        "input_path": str(path.resolve()),
        "checked_at_utc": utc_now(),
        "checks": checks,
        "available_external_tools": [
            record["name"] for record in tool_records if record.get("exists")
        ],
        "missing_external_tools": [
            record["name"] for record in tool_records if not record.get("exists")
        ],
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


def verify_idealized_summary(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    initial = payload.get("initial", {})
    final = payload.get("final", {})
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
    if schema_version == "gwm-next-boundary-cache/v1":
        return "boundary_cache_stub"
    if schema_version == "gwm-next-checkpoint-stub/v1":
        return "checkpoint_stub"
    if schema_version == "gwm-next-product-plan/v1":
        return "product_plan"
    if all(key in payload for key in ("case", "initial", "final")):
        return "idealized_summary"
    raise ValueError("Unable to detect input kind from JSON payload")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Verify prepared-case manifests, checkpoint stubs, and idealized summary products."
    )
    parser.add_argument("--input", required=True, help="JSON artifact path to verify")
    parser.add_argument(
        "--kind",
        default="auto",
        choices=[
            "auto",
            "prepared_case_manifest",
            "boundary_cache_stub",
            "checkpoint_stub",
            "product_plan",
            "idealized_summary",
        ],
    )
    parser.add_argument("--output-json", default=None, help="Optional explicit report path")
    args = parser.parse_args()

    input_path = Path(args.input)
    payload = load_json(input_path)
    kind = detect_kind(payload) if args.kind == "auto" else args.kind

    if kind == "prepared_case_manifest":
        report = verify_prepared_case_manifest(input_path, payload)
    elif kind == "boundary_cache_stub":
        report = verify_boundary_cache_stub(input_path, payload)
    elif kind == "checkpoint_stub":
        report = verify_checkpoint_stub(input_path, payload)
    elif kind == "product_plan":
        report = verify_product_plan(input_path, payload)
    elif kind == "idealized_summary":
        report = verify_idealized_summary(input_path, payload)
    else:  # pragma: no cover
        raise ValueError(f"Unsupported verification kind: {kind}")

    output_path = (
        Path(args.output_json)
        if args.output_json is not None
        else input_path.with_suffix(".verification.json")
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
