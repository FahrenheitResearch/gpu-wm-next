#include <filesystem>
#include <fstream>
#include <string>

#include "gwm/ingest/runtime_case.hpp"

#include "test_assert.hpp"

namespace {

std::filesystem::path write_text_file(const std::filesystem::path& path,
                                      const std::string& content) {
  std::ofstream out(path, std::ios::binary);
  out << content;
  return path;
}

}  // namespace

int main() {
  using namespace gwm::ingest;
  std::cerr << "runtime_case: start\n";

  const auto temp_dir = std::filesystem::temp_directory_path() / "gwm_runtime_case_test";
  std::filesystem::create_directories(temp_dir);
  std::cerr << "runtime_case: temp dir ready\n";

  const auto analysis_path = write_text_file(
      temp_dir / "analysis_state.json",
      R"({
  "schema_version": "gwm-next-analysis-state/v1",
  "source": "HRRR",
  "grid": {
    "nx": 2,
    "ny": 2,
    "nz": 2,
    "dx": 3000.0,
    "dy": 3000.0,
    "z_top": 12000.0
  },
  "cycle_time_utc": "2026-04-04T00:00:00Z",
  "valid_time_utc": "2026-04-04T00:00:00Z",
  "forecast_offset_seconds": 0,
  "metadata": {
    "status": "populated"
  },
  "atmosphere": {
    "u_wind": [1, 1, 1, 1, 2, 2, 2, 2],
    "v_wind": [0, 0, 0, 0, 1, 1, 1, 1],
    "w_wind": [0, 0, 0, 0, 0, 0, 0, 0],
    "air_temperature": [300, 300, 300, 300, 295, 295, 295, 295],
    "specific_humidity": [0.01, 0.01, 0.01, 0.01, 0.008, 0.008, 0.008, 0.008],
    "air_pressure": [90000, 90000, 90000, 90000, 80000, 80000, 80000, 80000],
    "geopotential_height": [100, 100, 100, 100, 1500, 1500, 1500, 1500]
  },
  "surface": {
    "surface_pressure": [95000, 95000, 95000, 95000],
    "air_temperature_2m": [299, 299, 299, 299],
    "specific_humidity_2m": [0.009, 0.009, 0.009, 0.009],
    "u_wind_10m": [8, 8, 8, 8],
    "v_wind_10m": [1, 1, 1, 1],
    "skin_temperature": [301, 301, 301, 301]
  },
  "static_surface": {
    "terrain_height": [100, 120, 140, 160],
    "land_mask": [1, 1, 0, 0],
    "land_use_index": [7, 7, 16, 16]
  }
}
)");
  std::cerr << "runtime_case: analysis fixture written\n";

  const auto boundary_path = write_text_file(
      temp_dir / "boundary_cache.json",
      R"({
  "schema_version": "gwm-next-boundary-cache/v1",
  "source": "HRRR",
  "grid": {
    "nx": 2,
    "ny": 2,
    "nz": 2,
    "dx": 3000.0,
    "dy": 3000.0,
    "z_top": 12000.0
  },
  "cycle_time_utc": "2026-04-04T00:00:00Z",
  "boundary_interval_seconds": 3600,
  "metadata": {
    "status": "populated"
  },
  "snapshots": [
    {
      "forecast_offset_seconds": 0,
      "valid_time_utc": "2026-04-04T00:00:00Z",
      "atmosphere": {
        "u_wind": [1, 1, 1, 1, 2, 2, 2, 2],
        "v_wind": [0, 0, 0, 0, 1, 1, 1, 1],
        "w_wind": [0, 0, 0, 0, 0, 0, 0, 0],
        "air_temperature": [300, 300, 300, 300, 295, 295, 295, 295],
        "specific_humidity": [0.01, 0.01, 0.01, 0.01, 0.008, 0.008, 0.008, 0.008],
        "air_pressure": [90000, 90000, 90000, 90000, 80000, 80000, 80000, 80000],
        "geopotential_height": [100, 100, 100, 100, 1500, 1500, 1500, 1500]
      },
      "surface": {
        "surface_pressure": [95000, 95000, 95000, 95000],
        "air_temperature_2m": [299, 299, 299, 299],
        "specific_humidity_2m": [0.009, 0.009, 0.009, 0.009],
        "u_wind_10m": [8, 8, 8, 8],
        "v_wind_10m": [1, 1, 1, 1],
        "skin_temperature": [301, 301, 301, 301]
      }
    },
    {
      "forecast_offset_seconds": 3600,
      "valid_time_utc": "2026-04-04T01:00:00Z",
      "atmosphere": {
        "u_wind": [3, 3, 3, 3, 4, 4, 4, 4],
        "v_wind": [1, 1, 1, 1, 2, 2, 2, 2],
        "w_wind": [0, 0, 0, 0, 0, 0, 0, 0],
        "air_temperature": [302, 302, 302, 302, 297, 297, 297, 297],
        "specific_humidity": [0.012, 0.012, 0.012, 0.012, 0.009, 0.009, 0.009, 0.009],
        "air_pressure": [90500, 90500, 90500, 90500, 80500, 80500, 80500, 80500],
        "geopotential_height": [110, 110, 110, 110, 1510, 1510, 1510, 1510]
      },
      "surface": {
        "surface_pressure": [95500, 95500, 95500, 95500],
        "air_temperature_2m": [300, 300, 300, 300],
        "specific_humidity_2m": [0.010, 0.010, 0.010, 0.010],
        "u_wind_10m": [9, 9, 9, 9],
        "v_wind_10m": [2, 2, 2, 2],
        "skin_temperature": [302, 302, 302, 302]
      }
    }
  ]
}
)");
  std::cerr << "runtime_case: boundary fixture written\n";

  std::cerr << "runtime_case: load analysis\n";
  const auto analysis = load_analysis_state_json(analysis_path.string());
  std::cerr << "runtime_case: analysis loaded\n";
  TEST_CHECK(analysis.source == SourceKind::HRRR);
  TEST_CHECK(analysis.grid.z_top == 12000.0);
  TEST_CHECK(analysis.atmosphere.has("u_wind"));

  std::cerr << "runtime_case: load cache\n";
  const auto cache = load_boundary_cache_json(boundary_path.string());
  std::cerr << "runtime_case: cache loaded\n";
  TEST_CHECK(cache.boundary_interval_seconds == 3600);
  TEST_CHECK(cache.snapshots.size() == 2);

  std::cerr << "runtime_case: load prepared case\n";
  const auto runtime_case =
      load_prepared_runtime_case(analysis_path.string(), boundary_path.string());
  std::cerr << "runtime_case: prepared case loaded\n";
  TEST_CHECK(runtime_case.analysis.grid.nx == 2);
  TEST_CHECK(runtime_case.boundary_cache.grid.nz == 2);

  std::cerr << "runtime_case: interpolate\n";
  const auto interpolated = interpolate_boundary_snapshot(cache, 1800);
  std::cerr << "runtime_case: interpolated\n";
  TEST_CHECK(interpolated.forecast_offset_seconds == 1800);
  TEST_NEAR(interpolated.atmosphere.values.at("u_wind")[0], 2.0f, 1.0e-6f);
  TEST_NEAR(interpolated.atmosphere.values.at("air_temperature")[0], 301.0f,
            1.0e-6f);
  TEST_NEAR(interpolated.atmosphere.values.at("air_pressure")[0], 90250.0f,
            1.0e-6f);
  TEST_NEAR(interpolated.surface.values.at("surface_pressure")[0], 95250.0f,
            1.0e-6f);

  std::filesystem::remove_all(temp_dir);
  std::cerr << "runtime_case: done\n";
  return 0;
}
