from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path


def parse_scalar(value: str):
    text = value.strip()
    if text.lower() in {"true", "false"}:
        return text.lower() == "true"
    try:
        if any(ch in text for ch in (".", "e", "E")):
            return float(text)
        return int(text)
    except ValueError:
        return text


def load_simple_yaml(path: Path) -> dict:
    data: dict[str, object] = {}
    active_section: str | None = None

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue

        indent = len(line) - len(line.lstrip(" "))
        stripped = line.strip()
        key, sep, value = stripped.partition(":")
        if not sep:
            raise ValueError(f"Invalid line in {path}: {raw_line}")

        if indent == 0:
            if value.strip():
                data[key] = parse_scalar(value)
                active_section = None
            else:
                data[key] = {}
                active_section = key
        elif indent == 2 and active_section is not None:
            assert isinstance(data[active_section], dict)
            data[active_section][key] = parse_scalar(value)
        else:
            raise ValueError(f"Unsupported YAML structure in {path}: {raw_line}")

    return data


def add_if_present(cmd: list[str], section: dict, key: str, flag: str | None = None):
    if key in section:
        cmd.extend([flag or f"--{key.replace('_', '-')}", str(section[key])])


def default_driver_path(repo_root: Path) -> Path:
    suffix = ".exe" if os.name == "nt" else ""
    return repo_root / "build-ninja" / f"gwm_idealized_driver{suffix}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run a gpu-wm-next idealized dry case from a simple case YAML."
    )
    parser.add_argument("case_file", type=Path)
    parser.add_argument(
        "--driver",
        type=Path,
        default=None,
        help="Path to gwm_idealized_driver executable",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=None,
        help="Optional explicit summary output path",
    )
    parser.add_argument("--steps", type=int, default=None)
    parser.add_argument("--dt", type=float, default=None)
    parser.add_argument("--fast-substeps", type=int, default=None)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    driver = args.driver or default_driver_path(repo_root)
    case_data = load_simple_yaml(args.case_file)

    grid = case_data.get("grid", {})
    terrain = case_data.get("terrain", {})
    init = case_data.get("initialization", {})
    if not isinstance(grid, dict) or not isinstance(terrain, dict) or not isinstance(init, dict):
        raise ValueError("Expected grid/terrain/initialization sections to be mappings")

    summary_json = args.summary_json or args.case_file.with_suffix(".summary.json")

    cmd = [str(driver), "--case", str(init.get("kind", "density_current"))]
    for key in ("nx", "ny", "nz", "dx", "dy", "z_top"):
        add_if_present(cmd, grid, key)

    if "kind" in terrain:
        cmd.extend(["--terrain-kind", str(terrain["kind"])])
        add_if_present(cmd, terrain, "center_x_fraction", "--terrain-center-x")
        add_if_present(cmd, terrain, "center_y_fraction", "--terrain-center-y")
        add_if_present(cmd, terrain, "half_width_x_fraction", "--terrain-half-width-x")
        add_if_present(cmd, terrain, "half_width_y_fraction", "--terrain-half-width-y")
        add_if_present(cmd, terrain, "height", "--terrain-height")

    for key in (
        "rho_background",
        "theta_background",
        "u_background",
        "v_background",
        "w_background",
        "theta_perturbation",
        "center_x_fraction",
        "center_y_fraction",
        "center_z_fraction",
        "radius_x_fraction",
        "radius_y_fraction",
        "radius_z_fraction",
        "rho_surface",
        "theta_ref",
        "density_scale_height",
    ):
        add_if_present(cmd, init, key)

    if args.steps is not None:
        cmd.extend(["--steps", str(args.steps)])
    if args.dt is not None:
        cmd.extend(["--dt", str(args.dt)])
    if args.fast_substeps is not None:
        cmd.extend(["--fast-substeps", str(args.fast_substeps)])

    cmd.extend(["--summary-json", str(summary_json)])

    print(f"Running idealized case: {args.case_file}")
    print(f"Driver: {driver}")
    print("Command:")
    print(" ".join(cmd))
    if os.name == "nt":
        wrapper = repo_root / "scripts" / "run_test_windows.cmd"
        subprocess.run(["cmd", "/c", str(wrapper), str(driver), *cmd[1:]], check=True)
    else:
        subprocess.run(cmd, check=True)
    print(f"Summary written to: {summary_json}")


if __name__ == "__main__":
    main()
