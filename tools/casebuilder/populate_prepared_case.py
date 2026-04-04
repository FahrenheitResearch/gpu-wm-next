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
KAPPA = 0.2854
DRY_CP = RD / KAPPA
LATENT_HEAT_VAPORIZATION = 2.5e6
WATER_VAPOR_GAS_CONSTANT = 461.5
EPSILON = RD / WATER_VAPOR_GAS_CONSTANT

LEGACY_SHORT_PRESSURE_LEVELS_HPA = [1000, 925, 850, 700, 500]
DEFAULT_DEEP_PRESSURE_LEVELS_HPA = [
    1000,
    975,
    950,
    925,
    900,
    875,
    850,
    825,
    800,
    775,
    750,
    725,
    700,
    675,
    650,
    625,
    600,
    575,
    550,
    525,
    500,
    475,
    450,
    425,
    400,
    375,
    350,
    325,
    300,
    275,
    250,
    225,
    200,
    175,
    150,
    125,
    100,
    75,
    50,
]


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


def terrain_weight_profile(eta: np.ndarray, terrain_taper_eta: float) -> np.ndarray:
    weights = np.empty_like(eta, dtype=np.float64)
    for index, eta_value in enumerate(eta):
        if eta_value <= 0.0:
            weights[index] = 1.0
        elif eta_value >= 1.0:
            weights[index] = 0.0
        elif terrain_taper_eta <= 0.0:
            s_value = float(np.clip(eta_value, 0.0, 1.0))
            weights[index] = 1.0 - s_value * s_value * (3.0 - 2.0 * s_value)
        elif eta_value <= terrain_taper_eta:
            weights[index] = 1.0
        else:
            s_value = float(
                np.clip(
                    (eta_value - terrain_taper_eta) / max(1.0e-9, 1.0 - terrain_taper_eta),
                    0.0,
                    1.0,
                )
            )
            weights[index] = 1.0 - s_value * s_value * (3.0 - 2.0 * s_value)
    return weights


def build_target_z_centers(
    nz: int,
    z_top: float,
    terrain_height: np.ndarray,
    terrain_taper_eta: float = 0.25,
) -> np.ndarray:
    eta_centers = (np.arange(nz, dtype=np.float64) + 0.5) / float(nz)
    flat_z_centers = eta_centers * float(z_top)
    terrain_weights = terrain_weight_profile(eta_centers, terrain_taper_eta)
    return (
        flat_z_centers[:, None, None]
        + terrain_weights[:, None, None] * terrain_height[None, :, :]
    ).astype(np.float32)


def interpolate_profile(
    source_heights: np.ndarray,
    source_values: np.ndarray,
    target_heights: np.ndarray,
    *,
    log_values: bool = False,
    clip_min: float | None = None,
) -> np.ndarray:
    order = np.argsort(source_heights)
    sorted_heights = np.asarray(source_heights[order], dtype=np.float64)
    sorted_values = np.asarray(source_values[order], dtype=np.float64)

    unique_heights, unique_indices = np.unique(sorted_heights, return_index=True)
    unique_values = sorted_values[unique_indices]

    if unique_heights.size == 1:
        result = np.full_like(target_heights, unique_values[0], dtype=np.float64)
    else:
        if log_values:
            transformed = np.log(np.maximum(unique_values, 1.0))
            result = np.exp(
                np.interp(
                    target_heights,
                    unique_heights,
                    transformed,
                    left=transformed[0],
                    right=transformed[-1],
                )
            )
        else:
            result = np.interp(
                target_heights,
                unique_heights,
                unique_values,
                left=unique_values[0],
                right=unique_values[-1],
            )

    if clip_min is not None:
        result = np.maximum(result, clip_min)
    return result.astype(np.float32)


def remap_atmosphere_to_model_levels(
    cropped_atmosphere: dict[str, np.ndarray],
    target_z_centers: np.ndarray,
) -> dict[str, np.ndarray]:
    nz_target, ny, nx = target_z_centers.shape
    source_heights = np.asarray(cropped_atmosphere["geopotential_height"], dtype=np.float32)
    remapped: dict[str, np.ndarray] = {
        "geopotential_height": np.array(target_z_centers, copy=True, dtype=np.float32)
    }
    linear_fields = ("u_wind", "v_wind", "w_wind", "air_temperature")
    nonnegative_fields = {"specific_humidity": 0.0}

    for field_name in linear_fields:
        remapped[field_name] = np.empty((nz_target, ny, nx), dtype=np.float32)
    for field_name in nonnegative_fields:
        remapped[field_name] = np.empty((nz_target, ny, nx), dtype=np.float32)
    remapped["air_pressure"] = np.empty((nz_target, ny, nx), dtype=np.float32)

    for j in range(ny):
        for i in range(nx):
            source_z_column = source_heights[:, j, i]
            target_z_column = target_z_centers[:, j, i]
            for field_name in linear_fields:
                remapped[field_name][:, j, i] = interpolate_profile(
                    source_z_column,
                    cropped_atmosphere[field_name][:, j, i],
                    target_z_column,
                )
            for field_name, clip_min in nonnegative_fields.items():
                remapped[field_name][:, j, i] = interpolate_profile(
                    source_z_column,
                    cropped_atmosphere[field_name][:, j, i],
                    target_z_column,
                    clip_min=clip_min,
                )
            remapped["air_pressure"][:, j, i] = interpolate_profile(
                source_z_column,
                cropped_atmosphere["air_pressure"][:, j, i],
                target_z_column,
                log_values=True,
                clip_min=1.0,
            )

    return remapped


def level_hpa_strings(levels_hpa: list[int]) -> str:
    return "|".join(str(level) for level in levels_hpa)


def normalize_pressure_levels(levels_hpa: list[int]) -> list[int]:
    return sorted({int(level) for level in levels_hpa}, reverse=True)


def choose_source_pressure_levels(
    requested_levels_hpa: list[int], target_nz: int
) -> list[int]:
    normalized = normalize_pressure_levels(requested_levels_hpa)
    if target_nz < 20 or len(normalized) >= 20:
        return normalized
    if normalized == LEGACY_SHORT_PRESSURE_LEVELS_HPA:
        return list(DEFAULT_DEEP_PRESSURE_LEVELS_HPA)
    return sorted(set(DEFAULT_DEEP_PRESSURE_LEVELS_HPA).union(normalized), reverse=True)


def saturation_vapor_pressure_pa(temperature_k: np.ndarray) -> np.ndarray:
    temp_c = temperature_k - 273.15
    es_hpa = 6.112 * np.exp((17.67 * temp_c) / np.maximum(temp_c + 243.5, 1.0e-6))
    return es_hpa * 100.0


def saturation_specific_humidity(
    temperature_k: np.ndarray, pressure_pa: np.ndarray
) -> np.ndarray:
    es_pa = saturation_vapor_pressure_pa(temperature_k)
    denom = np.maximum(pressure_pa - (1.0 - EPSILON) * es_pa, 1000.0)
    return np.clip(EPSILON * es_pa / denom, 0.0, 0.05)


def apply_warm_rain_saturation_adjustment(
    atmosphere: dict[str, np.ndarray],
) -> tuple[dict[str, np.ndarray], dict[str, float | int]]:
    temperature0 = np.asarray(atmosphere["air_temperature"], dtype=np.float64).copy()
    pressure = np.asarray(atmosphere["air_pressure"], dtype=np.float64)
    qv0 = np.asarray(atmosphere["specific_humidity"], dtype=np.float64).copy()
    qc0 = np.asarray(
        atmosphere.get("cloud_water_mixing_ratio", np.zeros_like(qv0)),
        dtype=np.float64,
    ).copy()
    qr = np.asarray(
        atmosphere.get("rain_water_mixing_ratio", np.zeros_like(qv0)),
        dtype=np.float64,
    ).copy()
    qt = np.clip(qv0, 0.0, None) + np.clip(qc0, 0.0, None)

    initial_rh = np.where(
        saturation_specific_humidity(temperature0, pressure) > 0.0,
        qv0 / np.maximum(saturation_specific_humidity(temperature0, pressure), 1.0e-12),
        0.0,
    )
    supersat_cells = int(np.count_nonzero(initial_rh > 1.0 + 1.0e-6))
    qv_eq = np.minimum(qt, saturation_specific_humidity(temperature0, pressure))
    for _ in range(32):
        temperature_eq = temperature0 + (LATENT_HEAT_VAPORIZATION / DRY_CP) * (
            qv0 - qv_eq
        )
        qsat_eq = saturation_specific_humidity(temperature_eq, pressure)
        new_qv_eq = np.minimum(qt, qsat_eq)
        if np.nanmax(np.abs(new_qv_eq - qv_eq)) <= 1.0e-9:
            qv_eq = new_qv_eq
            break
        qv_eq = new_qv_eq

    qv = np.clip(qv_eq, 0.0, None)
    qc = np.clip(qt - qv, 0.0, None)
    temperature = temperature0 + (LATENT_HEAT_VAPORIZATION / DRY_CP) * (qv0 - qv)
    total_condensed = float(np.sum(np.maximum(qc - qc0, 0.0)))

    adjusted = dict(atmosphere)
    adjusted["air_temperature"] = temperature.astype(np.float32)
    adjusted["specific_humidity"] = np.clip(qv, 0.0, None).astype(np.float32)
    adjusted["cloud_water_mixing_ratio"] = np.clip(qc, 0.0, None).astype(np.float32)
    adjusted["rain_water_mixing_ratio"] = np.clip(qr, 0.0, None).astype(np.float32)
    return adjusted, {
        "max_initial_relative_humidity_pct": float(np.nanmax(initial_rh) * 100.0),
        "supersaturated_cell_count": supersat_cells,
        "initial_condensed_mixing_ratio_sum": total_condensed,
    }


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
    requested_pressure_levels_hpa: list[int],
    valid_time_utc: str,
    forecast_offset_seconds: int,
) -> dict[str, Any]:
    pressure_levels_hpa = atmosphere_fields["pressure_levels_hpa"]
    grid = dict(prepared_case["grid"])
    grid["nx"] = crop_window.nx
    grid["ny"] = crop_window.ny
    target_nz = int(grid["nz"])
    terrain_height = crop_2d(surface_fields["terrain_height"], crop_window)
    cropped_atmosphere = {
        "u_wind": crop_3d(atmosphere_fields["u_wind"], crop_window),
        "v_wind": crop_3d(atmosphere_fields["v_wind"], crop_window),
        "w_wind": crop_3d(atmosphere_fields["w_wind"], crop_window),
        "air_temperature": crop_3d(atmosphere_fields["air_temperature"], crop_window),
        "specific_humidity": crop_3d(atmosphere_fields["specific_humidity"], crop_window),
        "air_pressure": crop_3d(atmosphere_fields["air_pressure"], crop_window),
        "geopotential_height": crop_3d(atmosphere_fields["geopotential_height"], crop_window),
    }
    source_top_height_m = float(np.nanmax(cropped_atmosphere["geopotential_height"]))
    target_z_centers = build_target_z_centers(
        target_nz,
        float(grid["z_top"]),
        terrain_height,
    )
    remapped_atmosphere = remap_atmosphere_to_model_levels(
        cropped_atmosphere,
        target_z_centers,
    )
    remapped_atmosphere, moist_adjustment = apply_warm_rain_saturation_adjustment(
        remapped_atmosphere
    )
    grid["ref_lat"] = float(np.mean(crop_2d(surface_fields["latitude"], crop_window)))
    grid["ref_lon"] = float(np.mean(crop_2d(surface_fields["longitude"], crop_window)))

    metadata = {
        "domain_name": prepared_case["domain_name"],
        "status": "populated",
        "requested_pressure_levels_hpa": normalize_pressure_levels(
            requested_pressure_levels_hpa
        ),
        "source_pressure_levels_hpa": pressure_levels_hpa,
        "source_top_height_m": source_top_height_m,
        "vertical_remap": {
            "target_nz": target_nz,
            "terrain_taper_eta": 0.25,
            "coordinate": "terrain-following-height-centers",
        },
        "warm_rain_initial_adjustment": moist_adjustment,
        "target_window": prepared_case.get("target_window", {}),
        "crop_window": {
            "i0": crop_window.i0,
            "i1": crop_window.i1,
            "j0": crop_window.j0,
            "j1": crop_window.j1,
        },
    }

    atmosphere = {
        "u_wind": flatten_3d(remapped_atmosphere["u_wind"]),
        "v_wind": flatten_3d(remapped_atmosphere["v_wind"]),
        "w_wind": flatten_3d(remapped_atmosphere["w_wind"]),
        "air_temperature": flatten_3d(remapped_atmosphere["air_temperature"]),
        "specific_humidity": flatten_3d(remapped_atmosphere["specific_humidity"]),
        "cloud_water_mixing_ratio": flatten_3d(
            remapped_atmosphere["cloud_water_mixing_ratio"]
        ),
        "rain_water_mixing_ratio": flatten_3d(
            remapped_atmosphere["rain_water_mixing_ratio"]
        ),
        "air_pressure": flatten_3d(remapped_atmosphere["air_pressure"]),
        "geopotential_height": flatten_3d(remapped_atmosphere["geopotential_height"]),
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
    requested_pressure_levels_hpa: list[int],
    source_pressure_levels_hpa: list[int],
    z_top: float,
) -> dict[str, Any]:
    grid = dict(prepared_case["grid"])
    grid["nx"] = crop_window.nx
    grid["ny"] = crop_window.ny
    grid["nz"] = int(grid["nz"])
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
            "requested_pressure_levels_hpa": normalize_pressure_levels(
                requested_pressure_levels_hpa
            ),
            "source_pressure_levels_hpa": normalize_pressure_levels(
                source_pressure_levels_hpa
            ),
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
    requested_pressure_levels_hpa = normalize_pressure_levels(
        [
            int(level)
            for level in target_window.get(
                "pressure_levels_hpa", LEGACY_SHORT_PRESSURE_LEVELS_HPA
            )
        ]
    )
    cycle_time_utc = prepared_case["times"]["cycle_time_utc"]
    grid = prepared_case["grid"]
    source_pressure_levels_hpa = choose_source_pressure_levels(
        requested_pressure_levels_hpa, int(grid["nz"])
    )

    surface_fields = extract_hrrr_surface_fields(cycle_time_utc, 0)
    atmosphere_fields = extract_hrrr_atmosphere_fields(
        cycle_time_utc, 0, source_pressure_levels_hpa
    )
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
        requested_pressure_levels_hpa,
        prepared_case["times"]["analysis_valid_time_utc"],
        0,
    )

    boundary_snapshots: list[dict[str, Any]] = []
    offsets = prepared_case["contracts"]["boundary_cache"]["forecast_offsets_seconds"]
    for offset in offsets:
        surface_snapshot = extract_hrrr_surface_fields(cycle_time_utc, int(offset))
        atmosphere_snapshot = extract_hrrr_atmosphere_fields(
            cycle_time_utc, int(offset), source_pressure_levels_hpa
        )
        analysis_like = build_analysis_payload(
            prepared_case,
            crop_window,
            surface_snapshot,
            atmosphere_snapshot,
            requested_pressure_levels_hpa,
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
        requested_pressure_levels_hpa,
        source_pressure_levels_hpa,
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
