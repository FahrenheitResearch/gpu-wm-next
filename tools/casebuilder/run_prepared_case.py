from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def plan_view_field_names(path: Path) -> list[str]:
    if not path.exists():
        return []
    payload = load_json(path)
    return [
        str(field.get("name", ""))
        for field in payload.get("fields", [])
        if field.get("name")
    ]


def runtime_summary_moisture(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    payload = load_json(path)
    final = payload.get("final", {})
    moisture = final.get("moisture", {})
    return moisture if isinstance(moisture, dict) else {}


def find_driver_binary(explicit: str | None) -> Path:
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit))

    for relative in (
        "build-ninja/gwm_prepared_case_driver.exe",
        "build/gwm_prepared_case_driver.exe",
        "build-vs/Release/gwm_prepared_case_driver.exe",
        "build-vs/Debug/gwm_prepared_case_driver.exe",
    ):
        candidates.append(REPO_ROOT / relative)

    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved.exists():
            return resolved
    raise FileNotFoundError(
        "Unable to find gwm_prepared_case_driver.exe; build the repo first or pass --driver-binary"
    )


def run_command(command: list[str], env: dict[str, str]) -> None:
    print("+", " ".join(command))
    subprocess.run(command, check=True, cwd=REPO_ROOT, env=env)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Populate and run a prepared-case manifest through the source-driven driver."
    )
    parser.add_argument(
        "--prepared-case",
        required=True,
        help="Path to prepared_case_manifest.json",
    )
    parser.add_argument(
        "--driver-binary",
        default=None,
        help="Optional explicit path to gwm_prepared_case_driver executable",
    )
    parser.add_argument("--populate", action="store_true")
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--dt", type=float, default=2.0)
    parser.add_argument("--fast-substeps", type=int, default=2)
    parser.add_argument("--plan-view-level", type=int, default=1)
    parser.add_argument("--render-maps", action="store_true")
    parser.add_argument("--map-output-dir", default=None)
    parser.add_argument(
        "--cuda-launch-blocking",
        type=int,
        choices=(0, 1),
        default=None,
        help="Optional debug override for CUDA_LAUNCH_BLOCKING. Default leaves the environment unchanged.",
    )
    parser.add_argument(
        "--map-fields",
        nargs="*",
        default=None,
        help="Optional subset of plan-view field names to render",
    )
    parser.add_argument("--skip-verify", action="store_true")
    args = parser.parse_args()

    prepared_case_path = Path(args.prepared_case).resolve()
    prepared_case = load_json(prepared_case_path)
    output_dir = prepared_case_path.parent

    env = os.environ.copy()
    env.setdefault("PYTHONIOENCODING", "utf-8")
    if args.cuda_launch_blocking is not None:
        env["CUDA_LAUNCH_BLOCKING"] = str(args.cuda_launch_blocking)

    if args.populate:
        run_command(
            [
                sys.executable,
                str(REPO_ROOT / "tools" / "casebuilder" / "populate_prepared_case.py"),
                "--prepared-case",
                str(prepared_case_path),
            ],
            env,
        )
        prepared_case = load_json(prepared_case_path)

    artifacts = prepared_case.get("artifacts", {})
    analysis_state = artifacts.get("analysis_state")
    boundary_cache = artifacts.get("boundary_cache")
    if not analysis_state or not boundary_cache:
        raise RuntimeError(
            "Prepared case is missing populated analysis_state/boundary_cache artifacts. "
            "Run with --populate or populate the manifest first."
        )

    summary_path = output_dir / "summary.json"
    plan_view_path = output_dir / "plan_view.json"
    driver_binary = find_driver_binary(args.driver_binary)

    run_command(
        [
            str(driver_binary),
            "--analysis-state",
            analysis_state,
            "--boundary-cache",
            boundary_cache,
            "--steps",
            str(args.steps),
            "--dt",
            str(args.dt),
            "--fast-substeps",
            str(args.fast_substeps),
            "--summary-json",
            str(summary_path),
            "--plan-view-json",
            str(plan_view_path),
            "--plan-view-level",
            str(args.plan_view_level),
        ],
        env,
    )

    emitted_fields = plan_view_field_names(plan_view_path)
    if args.map_fields:
        missing_fields = [
            field for field in args.map_fields if field not in emitted_fields
        ]
        if missing_fields:
            raise RuntimeError(
                "Requested map fields are not present in the emitted plan-view bundle: "
                f"{', '.join(missing_fields)}. Available fields: "
                f"{', '.join(emitted_fields) if emitted_fields else '<none>'}"
            )
    moist_fields = [
        field
        for field in emitted_fields
        if field
        in {
            "specific_humidity",
            "specific_humidity_2m",
            "cloud_water_mixing_ratio",
            "rain_water_mixing_ratio",
            "total_condensate",
            "column_cloud_water",
            "column_rain_water",
            "column_total_condensate",
            "column_rain_fraction",
            "accumulated_surface_precipitation",
            "mean_surface_precipitation_rate",
            "synthetic_reflectivity",
            "relative_humidity",
            "relative_humidity_2m",
            "dewpoint",
            "dewpoint_2m",
            "air_temperature",
            "air_temperature_2m",
            "air_pressure",
        }
    ]

    map_manifest_path: Path | None = None
    if args.render_maps:
        map_output_dir = (
            Path(args.map_output_dir).resolve()
            if args.map_output_dir
            else output_dir / "plan_view_maps"
        )
        render_command = [
            sys.executable,
            str(REPO_ROOT / "tools" / "verify" / "render_plan_view_maps.py"),
            "--input",
            str(plan_view_path),
            "--output-dir",
            str(map_output_dir),
        ]
        if args.map_fields:
            render_command.extend(["--fields", *args.map_fields])
        run_command(render_command, env)
        manifest_in_map_dir = map_output_dir / "map_manifest.json"
        if manifest_in_map_dir.exists():
            map_manifest_path = output_dir / "map_manifest.json"
            shutil.copyfile(manifest_in_map_dir, map_manifest_path)
        else:
            map_manifest_path = None

    if not args.skip_verify:
        run_command(
            [
                sys.executable,
                str(REPO_ROOT / "tools" / "verify" / "run_verification.py"),
                "--input",
                str(output_dir),
                "--kind",
                "source_run_bundle",
            ],
            env,
        )
        for artifact in (prepared_case_path, summary_path, plan_view_path):
            run_command(
                [
                    sys.executable,
                    str(REPO_ROOT / "tools" / "verify" / "run_verification.py"),
                    "--input",
                    str(artifact),
                ],
                env,
            )
        if map_manifest_path is not None:
            run_command(
                [
                    sys.executable,
                    str(REPO_ROOT / "tools" / "verify" / "run_verification.py"),
                    "--input",
                    str(map_manifest_path),
                ],
                env,
            )

    final_moisture = runtime_summary_moisture(summary_path)
    print(f"Prepared-case run complete: {prepared_case_path}")
    print(f"- summary: {summary_path}")
    print(f"- plan_view: {plan_view_path}")
    if emitted_fields:
        print(f"- plan_view_fields: {', '.join(emitted_fields)}")
    if moist_fields:
        print(f"- moist_fields: {', '.join(moist_fields)}")
    if final_moisture:
        print(
            "- final_moisture:"
            f" vapor={final_moisture.get('vapor_water_path_sum_kg_m2', 0.0)}"
            f" cloud={final_moisture.get('cloud_water_path_sum_kg_m2', 0.0)}"
            f" rain={final_moisture.get('rain_water_path_sum_kg_m2', 0.0)}"
            f" precip_sum_mm={final_moisture.get('accumulated_surface_precipitation_sum_mm', 0.0)}"
            f" precip_max_mm={final_moisture.get('max_surface_precipitation_mm', 0.0)}"
        )
    if map_manifest_path is not None:
        print(f"- map_manifest: {map_manifest_path}")


if __name__ == "__main__":
    main()
