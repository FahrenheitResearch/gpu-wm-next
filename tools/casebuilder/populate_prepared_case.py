from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

import numpy as np


RD = 287.05
GRAVITY = 9.81


@dataclass(frozen=True)
class CropWindow:
    i0: int
    i1: int
    j0: int
    j1: int

    @property
    def nx(self) -> int:
        return self.i1 - self.i0

    @property
    def ny(self) -> int:
        return self.j1 - self.j0


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def parse_utc(text: str) -> datetime:
    normalized = text.strip()
    if normalized.endswith("Z"):
        normalized = normalized[:-1] + "+00:00"
    return datetime.fromisoformat(normalized).astimezone(UTC)


def format_utc(value: datetime) -> str:
    return value.astimezone(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def utc_now() -> str:
    return format_utc(datetime.now(UTC))


def ensure_tool_paths(prepared_case: dict[str, Any]) -> None:
    toolchain = prepared_case.get("external_toolchain", {})
    manifest_path = toolchain.get("manifest_path")
    if manifest_path:
        manifest_dir = Path(manifest_path).resolve().parent
        if str(manifest_dir) not in sys.path:
            sys.path.insert(0, str(manifest_dir))

    for record in toolchain.get("tools", []):
        configured = record.get("configured_path", "")
        if not configured:
            continue
        repo_root = Path(configured)
        python_dir = repo_root / "python"
        if python_dir.exists():
            sys.path.insert(0, str(python_dir))


def choose_center_indices(
    latitude: np.ndarray,
    longitude: np.ndarray,
    center_lat: float | None,
    center_lon: float | None,
) -> tuple[int, int]:
    ny, nx = latitude.shape
    if center_lat is None or center_lon is None:
        return ny // 2, nx // 2

    lon_target = center_lon % 360.0
    dist2 = (latitude - center_lat) ** 2 + ((longitude % 360.0) - lon_target) ** 2
    flat_index = int(np.argmin(dist2))
    j_center, i_center = np.unravel_index(flat_index, dist2.shape)
    return int(j_center), int(i_center)


def compute_crop_window(
    latitude: np.ndarray,
    longitude: np.ndarray,
    target_nx: int,
    target_ny: int,
    center_lat: float | None,
    center_lon: float | None,
) -> CropWindow:
    ny, nx = latitude.shape
    if target_nx > nx or target_ny > ny:
        raise ValueError(
            f"Requested crop {target_nx}x{target_ny} exceeds source grid {nx}x{ny}"
        )

    j_center, i_center = choose_center_indices(latitude, longitude, center_lat, center_lon)
    i0 = max(0, i_center - target_nx // 2)
    j0 = max(0, j_center - target_ny // 2)
    i1 = min(nx, i0 + target_nx)
    j1 = min(ny, j0 + target_ny)
    i0 = i1 - target_nx
    j0 = j1 - target_ny
    return CropWindow(i0=i0, i1=i1, j0=j0, j1=j1)


def crop_2d(values: np.ndarray, window: CropWindow) -> np.ndarray:
    return np.asarray(values[window.j0 : window.j1, window.i0 : window.i1], dtype=np.float32)


def crop_3d(values: np.ndarray, window: CropWindow) -> np.ndarray:
    return np.asarray(values[:, window.j0 : window.j1, window.i0 : window.i1], dtype=np.float32)


def flatten_2d(values: np.ndarray) -> list[float]:
    return [float(value) for value in np.ravel(values, order="C")]


def flatten_3d(values: np.ndarray) -> list[float]:
    return [float(value) for value in np.ravel(values, order="C")]


def level_hpa_strings(levels_hpa: list[int]) -> str:
    return "|".join(str(level) for level in levels_hpa)


def open_source_dataset(
    model: str,
    product: str,
    cycle_time_utc: str,
    forecast_offset_seconds: int,
    search: str,
):
    from rusbie import Herbie

    fxx = int(round(forecast_offset_seconds / 3600.0))
    cycle_time_naive = parse_utc(cycle_time_utc).replace(tzinfo=None)
    handle = Herbie(cycle_time_naive, model=model, product=product, fxx=fxx, verbose=False)
    try:
        dataset = handle.xarray(search, remove_grib=False)
    except (FileNotFoundError, ValueError) as exc:
        if "No index file" not in str(exc) and not isinstance(exc, FileNotFoundError):
            raise
        handle.download(search, overwrite=False, errors="raise")
        handle = Herbie(cycle_time_naive, model=model, product=product, fxx=fxx, verbose=False)
        dataset = handle.xarray(search, remove_grib=False)
    if isinstance(dataset, list):
        if not dataset:
            raise RuntimeError(f"No datasets returned for search {search!r}")
        dataset = dataset[0]
    return dataset


def get_data_var(dataset, candidates: tuple[str, ...]):
    for name in candidates:
        if name in dataset.data_vars:
            return dataset[name]
    raise KeyError(f"None of the candidate variables {candidates!r} were present")


def pressure_velocity_to_w(omega_pa_s: np.ndarray, pressure_pa: np.ndarray, temperature_k: np.ndarray) -> np.ndarray:
    return -omega_pa_s * RD * temperature_k / np.maximum(pressure_pa * GRAVITY, 1.0)


def extract_first_data_var(
    cycle_time_utc: str,
    forecast_offset_seconds: int,
    search: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dataset = open_source_dataset(
        model="hrrr",
        product="sfc",
        cycle_time_utc=cycle_time_utc,
        forecast_offset_seconds=forecast_offset_seconds,
        search=search,
    )
    if not dataset.data_vars:
        raise RuntimeError(f"No data variables found for search {search!r}")
    var_name = next(iter(dataset.data_vars))
    values = np.asarray(dataset[var_name].values, dtype=np.float32)
    latitude = np.asarray(dataset["latitude"].values, dtype=np.float64)
    longitude = np.asarray(dataset["longitude"].values, dtype=np.float64)
    return values, latitude, longitude


def extract_hrrr_surface_fields(cycle_time_utc: str, forecast_offset_seconds: int) -> dict[str, np.ndarray]:
    surface_pressure, latitude, longitude = extract_first_data_var(
        cycle_time_utc, forecast_offset_seconds, ":PRES:surface:"
    )
    skin_temperature, _, _ = extract_first_data_var(
        cycle_time_utc, forecast_offset_seconds, ":TMP:surface:"
    )
    t2, _, _ = extract_first_data_var(
        cycle_time_utc, forecast_offset_seconds, ":TMP:2 m above ground:"
    )
    q2, _, _ = extract_first_data_var(
        cycle_time_utc, forecast_offset_seconds, ":SPFH:2 m above ground:"
    )
    u10, _, _ = extract_first_data_var(
        cycle_time_utc, forecast_offset_seconds, ":UGRD:10 m above ground:"
    )
    v10, _, _ = extract_first_data_var(
        cycle_time_utc, forecast_offset_seconds, ":VGRD:10 m above ground:"
    )
    terrain, _, _ = extract_first_data_var(
        cycle_time_utc, forecast_offset_seconds, ":HGT:surface:"
    )
    land_mask, _, _ = extract_first_data_var(
        cycle_time_utc, forecast_offset_seconds, ":LAND:surface:"
    )
    return {
        "latitude": latitude,
        "longitude": longitude,
        "surface_pressure": surface_pressure,
        "skin_temperature": skin_temperature,
        "air_temperature_2m": t2,
        "specific_humidity_2m": q2,
        "u_wind_10m": u10,
        "v_wind_10m": v10,
        "terrain_height": terrain,
        "land_mask": land_mask,
}


def extract_hrrr_prs_field(
    cycle_time_utc: str, forecast_offset_seconds: int, pressure_levels_hpa: list[int], field_name: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[int]]:
    level_pattern = level_hpa_strings(pressure_levels_hpa)
    dataset = open_source_dataset(
        model="hrrr",
        product="prs",
        cycle_time_utc=cycle_time_utc,
        forecast_offset_seconds=forecast_offset_seconds,
        search=f":{field_name}:({level_pattern}) mb:",
    )
    if not dataset.data_vars:
        raise RuntimeError(f"No pressure-level variables found for field {field_name!r}")
    var_name = next(iter(dataset.data_vars))
    values = np.asarray(dataset[var_name].values, dtype=np.float32)
    raw_levels = np.asarray(dataset["isobaricInhPa"].values, dtype=np.float32)
    pressure_hpa = raw_levels / 100.0 if float(np.nanmax(raw_levels)) > 2000.0 else raw_levels
    latitude = np.asarray(dataset["latitude"].values, dtype=np.float64)
    longitude = np.asarray(dataset["longitude"].values, dtype=np.float64)
    return values, latitude, longitude, [int(level) for level in pressure_hpa]


def extract_hrrr_atmosphere_fields(
    cycle_time_utc: str,
    forecast_offset_seconds: int,
    pressure_levels_hpa: list[int],
) -> dict[str, Any]:
    u, latitude, longitude, resolved_levels = extract_hrrr_prs_field(
        cycle_time_utc, forecast_offset_seconds, pressure_levels_hpa, "UGRD"
    )
    v, _, _, _ = extract_hrrr_prs_field(
        cycle_time_utc, forecast_offset_seconds, pressure_levels_hpa, "VGRD"
    )
    omega, _, _, _ = extract_hrrr_prs_field(
        cycle_time_utc, forecast_offset_seconds, pressure_levels_hpa, "VVEL"
    )
    temperature, _, _, _ = extract_hrrr_prs_field(
        cycle_time_utc, forecast_offset_seconds, pressure_levels_hpa, "TMP"
    )
    specific_humidity, _, _, _ = extract_hrrr_prs_field(
        cycle_time_utc, forecast_offset_seconds, pressure_levels_hpa, "SPFH"
    )
    geopotential_height, _, _, _ = extract_hrrr_prs_field(
        cycle_time_utc, forecast_offset_seconds, pressure_levels_hpa, "HGT"
    )
    pressure_hpa = np.asarray(resolved_levels, dtype=np.float32)
    pressure_pa = pressure_hpa[:, None, None] * 100.0
    w = pressure_velocity_to_w(omega, pressure_pa, temperature)

    return {
        "latitude": latitude,
        "longitude": longitude,
        "pressure_levels_hpa": [int(level) for level in pressure_hpa],
        "u_wind": u,
        "v_wind": v,
        "w_wind": w.astype(np.float32),
        "air_temperature": temperature,
        "specific_humidity": specific_humidity.astype(np.float32),
        "air_pressure": np.broadcast_to(
            pressure_pa.astype(np.float32), temperature.shape
        ).copy(),
        "geopotential_height": geopotential_height,
    }


def normalize_land_use(land_mask: np.ndarray) -> np.ndarray:
    return np.where(land_mask >= 0.5, 7.0, 16.0).astype(np.float32)


def build_analysis_payload(
    prepared_case: dict[str, Any],
    crop_window: CropWindow,
    surface_fields: dict[str, np.ndarray],
    atmosphere_fields: dict[str, Any],
    valid_time_utc: str,
    forecast_offset_seconds: int,
) -> dict[str, Any]:
    pressure_levels_hpa = atmosphere_fields["pressure_levels_hpa"]
    grid = dict(prepared_case["grid"])
    grid["nx"] = crop_window.nx
    grid["ny"] = crop_window.ny
    grid["nz"] = len(pressure_levels_hpa)
    cropped_geopotential_height = crop_3d(atmosphere_fields["geopotential_height"], crop_window)
    grid["z_top"] = max(
        float(grid.get("z_top", 0.0)),
        float(np.nanmax(cropped_geopotential_height)) * 1.15 + 250.0,
    )
    grid["ref_lat"] = float(np.mean(crop_2d(surface_fields["latitude"], crop_window)))
    grid["ref_lon"] = float(np.mean(crop_2d(surface_fields["longitude"], crop_window)))

    metadata = {
        "domain_name": prepared_case["domain_name"],
        "status": "populated",
        "pressure_levels_hpa": pressure_levels_hpa,
        "target_window": prepared_case.get("target_window", {}),
        "crop_window": {
            "i0": crop_window.i0,
            "i1": crop_window.i1,
            "j0": crop_window.j0,
            "j1": crop_window.j1,
        },
    }

    atmosphere = {
        "u_wind": flatten_3d(crop_3d(atmosphere_fields["u_wind"], crop_window)),
        "v_wind": flatten_3d(crop_3d(atmosphere_fields["v_wind"], crop_window)),
        "w_wind": flatten_3d(crop_3d(atmosphere_fields["w_wind"], crop_window)),
        "air_temperature": flatten_3d(crop_3d(atmosphere_fields["air_temperature"], crop_window)),
        "specific_humidity": flatten_3d(
            crop_3d(atmosphere_fields["specific_humidity"], crop_window)
        ),
        "air_pressure": flatten_3d(crop_3d(atmosphere_fields["air_pressure"], crop_window)),
        "geopotential_height": flatten_3d(cropped_geopotential_height),
    }
    surface = {
        "surface_pressure": flatten_2d(crop_2d(surface_fields["surface_pressure"], crop_window)),
        "air_temperature_2m": flatten_2d(crop_2d(surface_fields["air_temperature_2m"], crop_window)),
        "specific_humidity_2m": flatten_2d(crop_2d(surface_fields["specific_humidity_2m"], crop_window)),
        "u_wind_10m": flatten_2d(crop_2d(surface_fields["u_wind_10m"], crop_window)),
        "v_wind_10m": flatten_2d(crop_2d(surface_fields["v_wind_10m"], crop_window)),
        "skin_temperature": flatten_2d(crop_2d(surface_fields["skin_temperature"], crop_window)),
    }
    static_surface = {
        "terrain_height": flatten_2d(crop_2d(surface_fields["terrain_height"], crop_window)),
        "land_mask": flatten_2d(crop_2d(surface_fields["land_mask"], crop_window)),
        "land_use_index": flatten_2d(
            crop_2d(normalize_land_use(surface_fields["land_mask"]), crop_window)
        ),
    }

    return {
        "schema_version": "gwm-next-analysis-state/v1",
        "source": prepared_case["source"]["display_name"],
        "grid": grid,
        "cycle_time_utc": prepared_case["times"]["cycle_time_utc"],
        "valid_time_utc": valid_time_utc,
        "forecast_offset_seconds": forecast_offset_seconds,
        "metadata": metadata,
        "atmosphere": atmosphere,
        "surface": surface,
        "static_surface": static_surface,
    }


def build_boundary_payload(
    prepared_case: dict[str, Any],
    crop_window: CropWindow,
    snapshots: list[dict[str, Any]],
    pressure_levels_hpa: list[int],
    z_top: float,
) -> dict[str, Any]:
    grid = dict(prepared_case["grid"])
    grid["nx"] = crop_window.nx
    grid["ny"] = crop_window.ny
    grid["nz"] = len(pressure_levels_hpa)
    grid["z_top"] = z_top
    payload_snapshots: list[dict[str, Any]] = []
    for snapshot in snapshots:
        payload_snapshots.append(
            {
                "forecast_offset_seconds": snapshot["forecast_offset_seconds"],
                "valid_time_utc": snapshot["valid_time_utc"],
                "atmosphere": snapshot["atmosphere"],
                "surface": snapshot["surface"],
            }
        )

    return {
        "schema_version": "gwm-next-boundary-cache/v1",
        "source": prepared_case["source"]["display_name"],
        "grid": grid,
        "cycle_time_utc": prepared_case["times"]["cycle_time_utc"],
        "boundary_interval_seconds": prepared_case["source"]["boundary_interval_seconds"],
        "metadata": {
            "domain_name": prepared_case["domain_name"],
            "status": "populated",
            "pressure_levels_hpa": pressure_levels_hpa,
            "target_window": prepared_case.get("target_window", {}),
            "crop_window": {
                "i0": crop_window.i0,
                "i1": crop_window.i1,
                "j0": crop_window.j0,
                "j1": crop_window.j1,
            },
        },
        "snapshots": payload_snapshots,
    }


def populate_hrrr(prepared_case: dict[str, Any], output_dir: Path) -> dict[str, str]:
    target_window = prepared_case.get("target_window", {})
    pressure_levels_hpa = [int(level) for level in target_window.get("pressure_levels_hpa", [1000, 925, 850, 700, 500])]
    cycle_time_utc = prepared_case["times"]["cycle_time_utc"]
    grid = prepared_case["grid"]

    surface_fields = extract_hrrr_surface_fields(cycle_time_utc, 0)
    atmosphere_fields = extract_hrrr_atmosphere_fields(cycle_time_utc, 0, pressure_levels_hpa)
    crop_window = compute_crop_window(
        surface_fields["latitude"],
        surface_fields["longitude"],
        int(grid["nx"]),
        int(grid["ny"]),
        target_window.get("center_lat"),
        target_window.get("center_lon"),
    )

    analysis_payload = build_analysis_payload(
        prepared_case,
        crop_window,
        surface_fields,
        atmosphere_fields,
        prepared_case["times"]["analysis_valid_time_utc"],
        0,
    )

    boundary_snapshots: list[dict[str, Any]] = []
    offsets = prepared_case["contracts"]["boundary_cache"]["forecast_offsets_seconds"]
    for offset in offsets:
        surface_snapshot = extract_hrrr_surface_fields(cycle_time_utc, int(offset))
        atmosphere_snapshot = extract_hrrr_atmosphere_fields(cycle_time_utc, int(offset), pressure_levels_hpa)
        analysis_like = build_analysis_payload(
            prepared_case,
            crop_window,
            surface_snapshot,
            atmosphere_snapshot,
            valid_time_utc=format_utc(parse_utc(cycle_time_utc) + timedelta(seconds=int(offset))),
            forecast_offset_seconds=int(offset),
        )
        boundary_snapshots.append(
            {
                "forecast_offset_seconds": int(offset),
                "valid_time_utc": analysis_like["valid_time_utc"],
                "atmosphere": analysis_like["atmosphere"],
                "surface": analysis_like["surface"],
            }
        )

    boundary_payload = build_boundary_payload(
        prepared_case,
        crop_window,
        boundary_snapshots,
        analysis_payload["metadata"]["pressure_levels_hpa"],
        float(analysis_payload["grid"]["z_top"]),
    )

    analysis_path = output_dir / "analysis_state.json"
    boundary_path = output_dir / "boundary_cache.json"
    write_json(analysis_path, analysis_payload)
    write_json(boundary_path, boundary_payload)
    return {
        "analysis_state": str(analysis_path.resolve()),
        "boundary_cache": str(boundary_path.resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Populate a prepared-case manifest with real source data artifacts."
    )
    parser.add_argument(
        "--prepared-case",
        required=True,
        help="Path to prepared_case_manifest.json emitted by prepare_case.py",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional explicit directory for populated artifacts",
    )
    args = parser.parse_args()

    prepared_case_path = Path(args.prepared_case).resolve()
    prepared_case = load_json(prepared_case_path)
    ensure_tool_paths(prepared_case)

    source_key = prepared_case["source"]["key"]
    output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else prepared_case_path.parent
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    if source_key != "hrrr":
        raise SystemExit(
            f"populate_prepared_case currently supports source='hrrr' only, got {source_key!r}"
        )

    populated_paths = populate_hrrr(prepared_case, output_dir)

    prepared_case["artifacts"]["analysis_state"] = populated_paths["analysis_state"]
    prepared_case["artifacts"]["boundary_cache"] = populated_paths["boundary_cache"]
    prepared_case["external_toolchain"]["population_backend"] = {
        "downloader": "rusbie",
        "decoder": "cfrust",
    }
    prepared_case["populated_at_utc"] = utc_now()
    write_json(prepared_case_path, prepared_case)

    print(f"Populated prepared case: {prepared_case_path}")
    print(f"- analysis_state: {populated_paths['analysis_state']}")
    print(f"- boundary_cache: {populated_paths['boundary_cache']}")


if __name__ == "__main__":
    main()
