from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from datetime import UTC, datetime
from pathlib import Path


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Smoke-test source-run bundle verification.")
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep the generated bundle directory for inspection on failure.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    verifier = repo_root / "tools" / "verify" / "run_verification.py"
    cleanup = None
    if args.keep_temp:
        bundle_dir = Path(tempfile.mkdtemp(prefix="gwm_verify_source_run_"))
    else:
        temp_dir = tempfile.TemporaryDirectory(prefix="gwm_verify_source_run_")
        bundle_dir = Path(temp_dir.name)
        cleanup = temp_dir.cleanup
    try:
        image_dir = bundle_dir / "plan_view_maps"
        image_dir.mkdir(parents=True, exist_ok=True)
        image_path = image_dir / "rho_d.png"
        qv_image_path = image_dir / "specific_humidity.png"
        write_text(image_path, "placeholder image")
        write_text(qv_image_path, "placeholder image")

        prepared_case_manifest = {
            "schema_version": "gwm-next-prepared-case/v1",
            "prepared_at_utc": utc_now(),
            "populated_at_utc": utc_now(),
            "domain_name": "verification_smoke",
            "source": {
                "key": "hrrr",
                "display_name": "HRRR",
                "first_class": True,
                "boundary_interval_seconds": 3600,
            },
            "grid": {
                "nx": 2,
                "ny": 2,
                "nz": 2,
                "dx": 3000.0,
                "dy": 3000.0,
                "z_top": 12000.0,
            },
            "target_window": {
                "center_lat": 39.0,
                "center_lon": -97.0,
                "pressure_levels_hpa": [1000, 850],
            },
            "times": {
                "cycle_time_utc": "2026-04-04T00:00:00Z",
                "analysis_valid_time_utc": "2026-04-04T00:00:00Z",
                "forecast_hours": 1,
            },
            "contracts": {
                "analysis_state": {
                    "schema_version": "gwm-next-analysis-state/v1",
                    "required_field_groups": {
                        "atmosphere": [
                            "u_wind",
                            "v_wind",
                            "w_wind",
                            "air_temperature",
                            "specific_humidity",
                        ],
                        "surface": [
                            "surface_pressure",
                            "air_temperature_2m",
                            "specific_humidity_2m",
                            "u_wind_10m",
                            "v_wind_10m",
                            "skin_temperature",
                        ],
                        "static_surface": [
                            "terrain_height",
                            "land_mask",
                            "land_use_index",
                        ],
                    },
                },
                "boundary_cache": {
                    "schema_version": "gwm-next-boundary-cache/v1",
                    "interval_seconds": 3600,
                    "forecast_offsets_seconds": [0, 3600],
                },
                "checkpoint": {
                    "schema_version": "gwm-next-checkpoint-stub/v1",
                    "restart_family": "gwm-native-checkpoint",
                },
            },
            "external_toolchain": {
                "manifest_path": str(repo_root / "tools" / "casebuilder" / "toolchain_manifest.example.toml"),
                "manifest_exists": True,
                "policy": {
                    "allow_missing_optional_tools": True,
                    "require_runtime_core_separation": True,
                },
                "tools": [],
                "population_backend": {
                    "downloader": "rusbie",
                    "decoder": "cfrust",
                },
            },
            "artifacts": {
                "analysis_state": str(bundle_dir / "analysis_state.json"),
                "boundary_cache": str(bundle_dir / "boundary_cache.json"),
                "checkpoint_stub": str(bundle_dir / "checkpoint_stub.json"),
                "product_plan": str(bundle_dir / "product_plan.json"),
                "prepared_case_manifest": str(bundle_dir / "prepared_case_manifest.json"),
            },
        }

        analysis_state = {
            "schema_version": "gwm-next-analysis-state/v1",
            "source": "HRRR",
            "grid": {
                "nx": 2,
                "ny": 2,
                "nz": 2,
                "dx": 3000.0,
                "dy": 3000.0,
                "z_top": 12000.0,
            },
            "cycle_time_utc": "2026-04-04T00:00:00Z",
            "valid_time_utc": "2026-04-04T00:00:00Z",
            "forecast_offset_seconds": 0,
            "metadata": {
                "status": "populated",
                "target_window": {
                    "center_lat": 39.0,
                    "center_lon": -97.0,
                    "pressure_levels_hpa": [1000, 850],
                },
            },
            "atmosphere": {
                "u_wind": [1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0],
                "v_wind": [0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25],
                "w_wind": [0.0] * 8,
                "air_temperature": [300.0, 300.0, 300.0, 300.0, 295.0, 295.0, 295.0, 295.0],
                "specific_humidity": [0.01] * 8,
                "water_vapor_mixing_ratio": [0.0101] * 8,
                "air_pressure": [90000.0] * 4 + [80000.0] * 4,
                "geopotential_height": [100.0] * 4 + [1500.0] * 4,
            },
            "surface": {
                "surface_pressure": [95000.0, 95000.0, 95000.0, 95000.0],
                "air_temperature_2m": [299.0, 299.0, 299.0, 299.0],
                "specific_humidity_2m": [0.009] * 4,
                "u_wind_10m": [8.0] * 4,
                "v_wind_10m": [1.0] * 4,
                "skin_temperature": [301.0] * 4,
            },
            "static_surface": {
                "terrain_height": [100.0, 120.0, 140.0, 160.0],
                "land_mask": [1.0, 1.0, 0.0, 0.0],
                "land_use_index": [7.0, 7.0, 16.0, 16.0],
            },
        }

        boundary_cache = {
            "schema_version": "gwm-next-boundary-cache/v1",
            "source": "HRRR",
            "grid": analysis_state["grid"],
            "cycle_time_utc": "2026-04-04T00:00:00Z",
            "boundary_interval_seconds": 3600,
            "metadata": {
                "status": "populated",
                "target_window": prepared_case_manifest["target_window"],
            },
            "snapshots": [
                {
                    "forecast_offset_seconds": 0,
                    "valid_time_utc": "2026-04-04T00:00:00Z",
                    "atmosphere": analysis_state["atmosphere"],
                    "surface": analysis_state["surface"],
                },
                {
                    "forecast_offset_seconds": 3600,
                    "valid_time_utc": "2026-04-04T01:00:00Z",
                    "atmosphere": analysis_state["atmosphere"],
                    "surface": analysis_state["surface"],
                },
            ],
        }

        summary = {
            "case": "source_run_smoke",
            "steps": 1,
            "dt": 2.0,
            "initial": {
                "total_dry_mass": 8.0,
                "total_rho_theta_m": 2400.0,
                "tracers": {
                    "specific_humidity": {
                        "total_mass": 0.08,
                        "min": 0.01,
                        "max": 0.01,
                    }
                },
            },
            "final": {
                "total_dry_mass": 8.0,
                "total_rho_theta_m": 2400.0,
                "w_face": {"min": -0.1, "max": 0.1},
                "tracers": {
                    "specific_humidity": {
                        "total_mass": 0.08,
                        "min": 0.0095,
                        "max": 0.0105,
                    }
                },
            },
        }

        plan_view = {
            "schema_version": "gwm-next-plan-view/v1",
            "case": "source_run_smoke",
            "steps": 1,
            "dt": 2.0,
            "slice_k": 0,
            "grid": {
                "nx": 2,
                "ny": 2,
                "nz": 2,
                "dx": 3000.0,
                "dy": 3000.0,
            },
            "slice_mean_height_m": 50.0,
            "fields": [
                {
                    "name": "rho_d",
                    "units": "kg m^-3",
                    "location": "cell_center",
                    "nx": 2,
                    "ny": 2,
                    "storage": "row_major_yx",
                    "values": [1.0, 1.0, 1.0, 1.0],
                },
                {
                    "name": "specific_humidity",
                    "units": "kg kg^-1",
                    "location": "cell_center",
                    "nx": 2,
                    "ny": 2,
                    "storage": "row_major_yx",
                    "values": [0.010, 0.011, 0.009, 0.0105],
                },
                {
                    "name": "relative_humidity",
                    "units": "%",
                    "location": "cell_center",
                    "nx": 2,
                    "ny": 2,
                    "storage": "row_major_yx",
                    "values": [65.0, 70.0, 68.0, 72.0],
                },
                {
                    "name": "dewpoint",
                    "units": "K",
                    "location": "cell_center",
                    "nx": 2,
                    "ny": 2,
                    "storage": "row_major_yx",
                    "values": [293.0, 294.0, 292.5, 293.5],
                }
            ],
        }

        map_manifest = {
            "schema_version": "gwm-next-map-manifest/v1",
            "generated_at_utc": utc_now(),
            "input_plan_view_path": str(bundle_dir / "runtime_plan_view.json"),
            "case": "source_run_smoke",
            "slice_k": 0,
            "grid": {"nx": 2, "ny": 2},
            "images": [
                {
                    "field": "rho_d",
                    "units": "kg m^-3",
                    "path": str(image_path.resolve()),
                    "format": "png",
                },
                {
                    "field": "specific_humidity",
                    "units": "kg kg^-1",
                    "path": str(qv_image_path.resolve()),
                    "format": "png",
                }
            ],
        }

        checkpoint_stub = {
            "schema_version": "gwm-next-checkpoint-stub/v1",
            "restart_family": "gwm-native-checkpoint",
            "cold_start_contract": {
                "required_groups": ["atmosphere", "surface", "static_surface"],
            },
        }
        product_plan = {
            "schema_version": "gwm-next-product-plan/v1",
            "verification_modes": ["screen", "aloft", "products"],
            "planned_products": {
                "screen": ["surface_pressure"],
                "aloft": ["rho_d"],
                "products": ["verification_report_json"],
            },
        }

        write_json(bundle_dir / "prepared_case_manifest.json", prepared_case_manifest)
        write_json(bundle_dir / "analysis_state.json", analysis_state)
        write_json(bundle_dir / "boundary_cache.json", boundary_cache)
        write_json(bundle_dir / "checkpoint_stub.json", checkpoint_stub)
        write_json(bundle_dir / "product_plan.json", product_plan)
        write_json(bundle_dir / "runtime_summary.json", summary)
        write_json(bundle_dir / "runtime_plan_view.json", plan_view)
        write_json(bundle_dir / "plan_view_maps" / "map_manifest.json", map_manifest)

        subprocess.run(
            [
                sys.executable,
                str(verifier),
                "--input",
                str(bundle_dir),
                "--kind",
                "source_run_bundle",
            ],
            cwd=repo_root,
            check=True,
        )

        subprocess.run(
            [
                sys.executable,
                str(verifier),
                "--input",
                str(bundle_dir / "prepared_case_manifest.json"),
            ],
            cwd=repo_root,
            check=True,
        )

        subprocess.run(
            [
                sys.executable,
                str(verifier),
                "--input",
                str(bundle_dir / "analysis_state.json"),
            ],
            cwd=repo_root,
            check=True,
        )

        subprocess.run(
            [
                sys.executable,
                str(verifier),
                "--input",
                str(bundle_dir / "boundary_cache.json"),
            ],
            cwd=repo_root,
            check=True,
        )

        subprocess.run(
            [
                sys.executable,
                str(verifier),
                "--input",
                str(bundle_dir / "runtime_summary.json"),
            ],
            cwd=repo_root,
            check=True,
        )

        subprocess.run(
            [
                sys.executable,
                str(verifier),
                "--input",
                str(bundle_dir / "runtime_plan_view.json"),
            ],
            cwd=repo_root,
            check=True,
        )

        subprocess.run(
            [
                sys.executable,
                str(verifier),
                "--input",
                str(bundle_dir / "plan_view_maps" / "map_manifest.json"),
            ],
            cwd=repo_root,
            check=True,
        )

        print("source-run bundle verification smoke passed")
    finally:
        if args.keep_temp:
            print(f"kept temp bundle at: {bundle_dir}")
        elif cleanup is not None:
            cleanup()


if __name__ == "__main__":
    main()
