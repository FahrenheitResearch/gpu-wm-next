#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/core/dry_thermo.hpp"
#include "gwm/dycore/passive_tracer.hpp"
#include "gwm/io/plan_view_output.hpp"

#include "test_assert.hpp"

namespace {

const gwm::io::PlanViewField& require_field(const gwm::io::PlanViewBundle& bundle,
                                            const std::string& name) {
  const auto it = std::find_if(
      bundle.fields.begin(), bundle.fields.end(),
      [&](const auto& field) { return field.name == name; });
  TEST_CHECK(it != bundle.fields.end());
  return *it;
}

int count_occurrences(const std::string& text, const std::string& needle) {
  if (needle.empty()) {
    return 0;
  }

  int count = 0;
  std::size_t pos = 0;
  while ((pos = text.find(needle, pos)) != std::string::npos) {
    ++count;
    pos += needle.size();
  }
  return count;
}

std::filesystem::path write_text_file(const std::filesystem::path& path,
                                      const std::string& content) {
  std::ofstream out(path, std::ios::binary);
  out << content;
  return path;
}

}  // namespace

int main() {
  try {
    using namespace gwm;

    const auto layout =
        comm::VirtualRankLayout::build(6, 4, 3, 1, 2, 1, true, true);
    auto states = dycore::make_constant_dry_state(layout, 1.2f, 300.0f, 10.0f,
                                                  -4.0f, 0.5f,
                                                  "plan_view_test");
    const auto metrics = domain::GridMetrics::make_hybrid_height(
        6, 4, 3, 1000.0, 1000.0, 9000.0, std::vector<real>(24, 0.0f));

    const auto bundle = io::extract_dry_plan_view(states, layout, metrics,
                                                  "constant", 4, 2.0f, 1);
    TEST_CHECK(bundle.schema_version == "gwm-next-plan-view/v1");
    TEST_CHECK(bundle.case_kind == "constant");
    TEST_CHECK(bundle.steps == 4);
    TEST_NEAR(bundle.dt, 2.0f, 1.0e-6f);
    TEST_CHECK(bundle.slice_k == 1);
    TEST_CHECK(bundle.nx == 6);
    TEST_CHECK(bundle.ny == 4);
    TEST_NEAR(bundle.dx, 1000.0, 1.0e-9);
    TEST_NEAR(bundle.dy, 1000.0, 1.0e-9);
    TEST_NEAR(bundle.slice_mean_height_m, 4500.0, 1.0e-6);
    TEST_CHECK(bundle.fields.size() == 10);

    std::set<std::string> unique_names;
    for (const auto& field : bundle.fields) {
      TEST_CHECK(unique_names.insert(field.name).second);
      TEST_CHECK(field.location == "cell_center");
      TEST_CHECK(!field.units.empty());
      TEST_CHECK(field.nx == bundle.nx);
      TEST_CHECK(field.ny == bundle.ny);
      TEST_CHECK(static_cast<int>(field.values.size()) == bundle.nx * bundle.ny);
    }

    const auto& terrain = require_field(bundle, "terrain_height");
    const auto& z_center = require_field(bundle, "z_center");
    const auto& rho = require_field(bundle, "rho_d");
    const auto& theta = require_field(bundle, "theta_m");
    const auto& pressure = require_field(bundle, "air_pressure");
    const auto& temperature = require_field(bundle, "air_temperature");
    const auto& u = require_field(bundle, "u_velocity");
    const auto& v = require_field(bundle, "v_velocity");
    const auto& wind = require_field(bundle, "wind_speed");
    const auto& w = require_field(bundle, "w_velocity");

    TEST_CHECK(static_cast<int>(terrain.values.size()) == bundle.nx * bundle.ny);
    TEST_CHECK(static_cast<int>(rho.values.size()) == bundle.nx * bundle.ny);
    TEST_NEAR(rho.values.front(), 1.2, 1.0e-6);
    TEST_NEAR(theta.values.front(), 300.0, 1.0e-6);
    TEST_NEAR(u.values.front(), 10.0, 1.0e-6);
    TEST_NEAR(v.values.front(), -4.0, 1.0e-6);
    TEST_NEAR(w.values.front(), 0.5, 1.0e-6);
    TEST_NEAR(wind.values.front(), std::sqrt(116.0), 1.0e-6);
    TEST_NEAR(terrain.values.front(), 0.0, 1.0e-6);
    TEST_NEAR(z_center.values.front(), 4500.0, 1.0e-6);
    const auto expected_pressure =
        core::dry_pressure_from_rho_theta_m(1.2f, 1.2f * 300.0f);
    const auto expected_temperature =
        300.0f *
        std::pow(expected_pressure / core::kReferencePressure, core::kKappa);
    TEST_NEAR(pressure.values.front(), expected_pressure, 1.0e-4);
    TEST_NEAR(temperature.values.front(), expected_temperature, 1.0e-4);

    const std::array<const char*, 10> expected_fields = {
        "terrain_height", "z_center",       "rho_d",          "theta_m",
        "air_pressure",   "air_temperature","u_velocity",     "v_velocity",
        "wind_speed",     "w_velocity"};
    for (const auto* expected : expected_fields) {
      TEST_CHECK(unique_names.count(expected) == 1);
    }

    const auto json = io::plan_view_bundle_to_json(bundle);
    TEST_CHECK(json.find("\"schema_version\": \"gwm-next-plan-view/v1\"") !=
               std::string::npos);
    TEST_CHECK(json.find("\"case\": \"constant\"") != std::string::npos);
    TEST_CHECK(json.find("\"slice_k\": 1") != std::string::npos);
    TEST_CHECK(json.find("\"grid\": {") != std::string::npos);
    TEST_CHECK(json.find("\"slice_mean_height_m\": 4500") != std::string::npos);
    TEST_CHECK(json.find("\"name\": \"air_pressure\"") != std::string::npos);
    TEST_CHECK(json.find("\"name\": \"air_temperature\"") != std::string::npos);
    TEST_CHECK(json.find("\"name\": \"theta_m\"") != std::string::npos);
    TEST_CHECK(json.find("\"storage\": \"row_major_yx\"") != std::string::npos);
    TEST_CHECK(count_occurrences(json, "\"storage\": \"row_major_yx\"") == 10);
    TEST_CHECK(count_occurrences(json, "\"location\": \"cell_center\"") == 10);

    std::vector<real> qv_global(static_cast<std::size_t>(metrics.nx) *
                                    static_cast<std::size_t>(metrics.ny) *
                                    static_cast<std::size_t>(metrics.nz),
                                0.004f);
    const auto runtime_tracers =
        dycore::make_specific_humidity_tracers_from_global_field(
            states, layout, metrics.nx, metrics.ny, qv_global, "runtime_qv");
    const auto runtime_bundle = io::extract_runtime_plan_view(
        states, runtime_tracers, layout, metrics, "constant_runtime", 4, 2.0f, 1);
    TEST_CHECK(runtime_bundle.fields.size() == 13);
    const auto& runtime_q = require_field(runtime_bundle, "specific_humidity");
    const auto& runtime_rh = require_field(runtime_bundle, "relative_humidity");
    const auto& runtime_dewpoint = require_field(runtime_bundle, "dewpoint");
    TEST_NEAR(runtime_q.values.front(), 0.004, 1.0e-6);
    TEST_CHECK(runtime_rh.values.front() > 0.0);
    TEST_CHECK(runtime_rh.values.front() <= 100.0);
    TEST_CHECK(runtime_dewpoint.values.front() < temperature.values.front());

    auto enriched = bundle;
    const auto temp_dir =
        std::filesystem::temp_directory_path() / "gwm_plan_view_output_test";
    std::filesystem::create_directories(temp_dir);
    const auto analysis_path = write_text_file(
        temp_dir / "analysis_state.json",
        R"({
  "schema_version": "gwm-next-analysis-state/v1",
  "source": "HRRR",
  "grid": {
    "nx": 6,
    "ny": 4,
    "nz": 3,
    "dx": 1000.0,
    "dy": 1000.0,
    "z_top": 9000.0
  },
  "cycle_time_utc": "2026-04-04T00:00:00Z",
  "valid_time_utc": "2026-04-04T00:00:00Z",
  "forecast_offset_seconds": 0,
  "metadata": {
    "status": "populated"
  },
  "atmosphere": {
    "u_wind": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "v_wind": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "w_wind": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "air_temperature": [300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300],
    "specific_humidity": [0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014],
    "air_pressure": [100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 90000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000, 80000],
    "geopotential_height": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000]
  },
  "surface": {
    "surface_pressure": [100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000],
    "air_temperature_2m": [300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300],
    "specific_humidity_2m": [0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004],
    "u_wind_10m": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "v_wind_10m": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "skin_temperature": [300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300]
  },
  "static_surface": {
    "terrain_height": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "land_mask": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    "land_use_index": [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
  }
}
)");

    io::enrich_plan_view_bundle_from_prepared_case(enriched, analysis_path.string());
    TEST_CHECK(enriched.fields.size() == 18);
    const auto& q = require_field(enriched, "specific_humidity");
    const auto& rh = require_field(enriched, "relative_humidity");
    const auto& dewpoint = require_field(enriched, "dewpoint");
    TEST_CHECK(require_field(enriched, "surface_pressure").values.size() ==
               static_cast<std::size_t>(bundle.nx * bundle.ny));
    TEST_CHECK(require_field(enriched, "air_temperature_2m").values.size() ==
               static_cast<std::size_t>(bundle.nx * bundle.ny));
    TEST_CHECK(require_field(enriched, "specific_humidity_2m").values.size() ==
               static_cast<std::size_t>(bundle.nx * bundle.ny));
    TEST_CHECK(require_field(enriched, "relative_humidity_2m").values.size() ==
               static_cast<std::size_t>(bundle.nx * bundle.ny));
    TEST_CHECK(require_field(enriched, "dewpoint_2m").values.size() ==
               static_cast<std::size_t>(bundle.nx * bundle.ny));
    TEST_NEAR(q.values.front(), 0.009, 1.0e-6);
    TEST_CHECK(rh.values.front() > 0.0);
    TEST_CHECK(rh.values.front() <= 100.0);
    TEST_CHECK(dewpoint.values.front() < temperature.values.front());

    const auto output_path = temp_dir / "plan_view.json";
    io::write_plan_view_bundle_json(bundle, output_path.string());
    std::ifstream in(output_path, std::ios::binary);
    const std::string enriched_json((std::istreambuf_iterator<char>(in)),
                                    std::istreambuf_iterator<char>());
    TEST_CHECK(enriched_json.find("\"name\": \"specific_humidity\"") !=
               std::string::npos);
    TEST_CHECK(enriched_json.find("\"name\": \"relative_humidity\"") !=
               std::string::npos);
    TEST_CHECK(enriched_json.find("\"name\": \"dewpoint\"") !=
               std::string::npos);
    TEST_CHECK(enriched_json.find("\"name\": \"specific_humidity_2m\"") !=
               std::string::npos);
    in.close();

    std::filesystem::remove_all(temp_dir);
    return 0;
  } catch (const std::exception& ex) {
    test_fail(std::string("exception: ") + ex.what());
  }
}
