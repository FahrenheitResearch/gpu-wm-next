from __future__ import annotations

import math
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]

import sys

sys.path.insert(0, str(REPO_ROOT / "tools" / "casebuilder"))

import populate_prepared_case as pop


def assert_close(actual: float, expected: float, tol: float = 1.0e-6) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"expected {expected}, got {actual}")


def test_choose_source_pressure_levels() -> None:
    expanded = pop.choose_source_pressure_levels(
        [1000, 925, 850, 700, 500], target_nz=50
    )
    if expanded == [1000, 925, 850, 700, 500]:
        raise AssertionError("large-nz source levels were not expanded")
    if 100 not in expanded or 50 not in expanded:
        raise AssertionError("expanded source levels do not reach the upper troposphere")
    if expanded[0] != 1000 or expanded[-1] != 50:
        raise AssertionError("expanded source levels are not ordered as expected")

    preserved = pop.choose_source_pressure_levels([1000, 850], target_nz=2)
    if preserved != [1000, 850]:
        raise AssertionError("small-nz case should preserve the requested level list")


def test_saturation_adjustment() -> None:
    temperature = np.full((2, 1, 1), 280.0, dtype=np.float32)
    pressure = np.full((2, 1, 1), 90000.0, dtype=np.float32)
    qv = np.full((2, 1, 1), 0.02, dtype=np.float32)
    adjusted, stats = pop.apply_warm_rain_saturation_adjustment(
        {
            "air_temperature": temperature,
            "air_pressure": pressure,
            "specific_humidity": qv,
        }
    )

    final_qv = adjusted["specific_humidity"]
    final_qc = adjusted["cloud_water_mixing_ratio"]
    final_t = adjusted["air_temperature"]
    final_qsat = pop.saturation_specific_humidity(final_t, pressure)

    if not np.any(final_qc > 0.0):
        raise AssertionError("supersaturated source column did not create cloud water")
    if not np.all(final_qv <= final_qsat + 1.0e-5):
        raise AssertionError("saturation adjustment left the column supersaturated")
    if not np.all(final_t >= temperature):
        raise AssertionError("saturation adjustment did not warm the supersaturated cells")
    if stats["supersaturated_cell_count"] <= 0:
        raise AssertionError("saturation-adjustment metadata did not record supersaturation")


def test_build_analysis_preserves_manifest_z_top() -> None:
    prepared_case = {
        "domain_name": "unit_test_case",
        "source": {"display_name": "HRRR"},
        "grid": {
            "nx": 2,
            "ny": 2,
            "nz": 4,
            "dx": 3000.0,
            "dy": 3000.0,
            "z_top": 12000.0,
        },
        "times": {"cycle_time_utc": "2026-04-04T00:00:00Z"},
        "target_window": {
            "center_lat": 39.0,
            "center_lon": -97.0,
            "pressure_levels_hpa": [1000, 850],
        },
    }
    window = pop.CropWindow(i0=0, i1=2, j0=0, j1=2)
    surface_fields = {
        "latitude": np.array([[39.0, 39.0], [39.1, 39.1]], dtype=np.float64),
        "longitude": np.array([[-97.0, -96.9], [-97.0, -96.9]], dtype=np.float64),
        "surface_pressure": np.full((2, 2), 95000.0, dtype=np.float32),
        "skin_temperature": np.full((2, 2), 301.0, dtype=np.float32),
        "air_temperature_2m": np.full((2, 2), 299.0, dtype=np.float32),
        "specific_humidity_2m": np.full((2, 2), 0.009, dtype=np.float32),
        "u_wind_10m": np.full((2, 2), 8.0, dtype=np.float32),
        "v_wind_10m": np.full((2, 2), 1.0, dtype=np.float32),
        "terrain_height": np.array([[100.0, 120.0], [140.0, 160.0]], dtype=np.float32),
        "land_mask": np.array([[1.0, 1.0], [0.0, 0.0]], dtype=np.float32),
    }
    shape = (2, 2, 2)
    atmosphere_fields = {
        "pressure_levels_hpa": [1000, 850],
        "u_wind": np.full(shape, 5.0, dtype=np.float32),
        "v_wind": np.full(shape, 1.0, dtype=np.float32),
        "w_wind": np.zeros(shape, dtype=np.float32),
        "air_temperature": np.array(
            [
                [[290.0, 290.0], [290.0, 290.0]],
                [[275.0, 275.0], [275.0, 275.0]],
            ],
            dtype=np.float32,
        ),
        "specific_humidity": np.array(
            [
                [[0.016, 0.016], [0.016, 0.016]],
                [[0.004, 0.004], [0.004, 0.004]],
            ],
            dtype=np.float32,
        ),
        "air_pressure": np.array(
            [
                [[90000.0, 90000.0], [90000.0, 90000.0]],
                [[70000.0, 70000.0], [70000.0, 70000.0]],
            ],
            dtype=np.float32,
        ),
        "geopotential_height": np.array(
            [
                [[500.0, 500.0], [500.0, 500.0]],
                [[14000.0, 14000.0], [14000.0, 14000.0]],
            ],
            dtype=np.float32,
        ),
    }

    payload = pop.build_analysis_payload(
        prepared_case,
        window,
        surface_fields,
        atmosphere_fields,
        requested_pressure_levels_hpa=[1000, 850],
        valid_time_utc="2026-04-04T00:00:00Z",
        forecast_offset_seconds=0,
    )

    assert_close(payload["grid"]["z_top"], 12000.0)
    if payload["metadata"]["requested_pressure_levels_hpa"] != [1000, 850]:
        raise AssertionError("requested pressure-level metadata was not preserved")
    if payload["metadata"]["source_pressure_levels_hpa"] != [1000, 850]:
        raise AssertionError("source pressure-level metadata was not recorded")
    if payload["metadata"]["source_top_height_m"] <= payload["grid"]["z_top"]:
        raise AssertionError("source-top metadata should reflect the taller input column")
    atmosphere = payload["atmosphere"]
    if "cloud_water_mixing_ratio" not in atmosphere or "rain_water_mixing_ratio" not in atmosphere:
        raise AssertionError("warm-rain tracer fields were not emitted into analysis payload")


def main() -> None:
    test_choose_source_pressure_levels()
    test_saturation_adjustment()
    test_build_analysis_preserves_manifest_z_top()


if __name__ == "__main__":
    main()
