#include <cstdlib>
#include <exception>
#include <string>
#include <vector>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/domain/grid_metrics.hpp"
#include "gwm/core/runtime_summary.hpp"
#include "gwm/state/tracer_registry.hpp"

#include "test_assert.hpp"

namespace {

const gwm::core::TracerSummary& require_tracer(
    const gwm::core::RuntimeStateSummary& summary, const std::string& name) {
  for (const auto& tracer : summary.tracers) {
    if (tracer.name == name) {
      return tracer;
    }
  }
  test_fail("missing tracer in runtime summary: " + name);
  std::abort();
}

}  // namespace

int main() {
  try {
    using namespace gwm;

    std::vector<dycore::DryState> states;
    states.emplace_back(2, 1, 1, 1, "runtime_summary");
    states.front().fill_constant(1.0f, 300.0f, 0.0f, 0.0f, 0.0f);
    const auto layout = comm::VirtualRankLayout::build(2, 1, 1, 1, 1, 1);
    const auto metrics = domain::GridMetrics::make_hybrid_height(
        2, 1, 1, 1000.0, 1000.0, 1000.0,
        std::vector<real>(static_cast<std::size_t>(2), 0.0f));

    std::vector<state::TracerState> tracers;
    tracers.emplace_back(state::make_warm_rain_registry(), 2, 1, 1, 1,
                         "runtime_summary");
    tracers.front().fill_zero();
    tracers.front().mass(state::kSpecificHumidityTracerName)(0, 0, 0) = 0.010f;
    tracers.front().mass(state::kSpecificHumidityTracerName)(1, 0, 0) = 0.012f;
    tracers.front().mass(state::kCloudWaterTracerName)(0, 0, 0) = 0.001f;
    tracers.front().mass(state::kRainWaterTracerName)(1, 0, 0) = 0.002f;
    std::vector<physics::WarmRainSurfaceAccumulation> accumulations(1);
    accumulations.front().reset(2, 1);
    accumulations.front().liquid_precipitation_kg_m2(1, 0, 0) = 0.75f;
    const double dz0 =
        1.0 / static_cast<double>(metrics.inv_dz_cell(0, 0, 0));
    const double dz1 =
        1.0 / static_cast<double>(metrics.inv_dz_cell(1, 0, 0));
    const double expected_qv_path = 0.010 * dz0 + 0.012 * dz1;
    const double expected_qc_path = 0.001 * dz0;
    const double expected_qr_path = 0.002 * dz1;

    const auto summary =
        core::summarize_runtime_state(states, tracers, layout, metrics,
                                      &accumulations);
    const auto& qv = require_tracer(summary, state::kSpecificHumidityTracerName);
    const auto& qc = require_tracer(summary, state::kCloudWaterTracerName);
    const auto& qr = require_tracer(summary, state::kRainWaterTracerName);
    TEST_NEAR(summary.moisture.vapor_water_path_sum_kg_m2, expected_qv_path,
              1.0e-5);
    TEST_NEAR(summary.moisture.cloud_water_path_sum_kg_m2, expected_qc_path,
              1.0e-5);
    TEST_NEAR(summary.moisture.rain_water_path_sum_kg_m2, expected_qr_path,
              1.0e-5);
    TEST_NEAR(summary.moisture.condensed_water_path_sum_kg_m2,
              expected_qc_path + expected_qr_path, 1.0e-5);
    TEST_NEAR(summary.moisture.total_water_path_sum_kg_m2,
              expected_qv_path + expected_qc_path + expected_qr_path, 1.0e-5);
    TEST_NEAR(summary.moisture.accumulated_surface_precipitation_sum_mm, 0.75,
              1.0e-6);
    TEST_NEAR(summary.moisture.mean_surface_precipitation_mm, 0.375, 1.0e-6);
    TEST_NEAR(summary.moisture.max_surface_precipitation_mm, 0.75, 1.0e-6);
    TEST_CHECK(summary.tracers.size() == 3);
    TEST_NEAR(qv.mixing_ratio.min, 0.010, 1.0e-6);
    TEST_NEAR(qv.mixing_ratio.max, 0.012, 1.0e-6);
    TEST_NEAR(qc.mixing_ratio.min, 0.0, 1.0e-6);
    TEST_NEAR(qc.mixing_ratio.max, 0.001, 1.0e-6);
    TEST_NEAR(qr.mixing_ratio.min, 0.0, 1.0e-6);
    TEST_NEAR(qr.mixing_ratio.max, 0.002, 1.0e-6);
    TEST_CHECK(qv.positive);
    TEST_CHECK(qc.positive);
    TEST_CHECK(qr.positive);

    const auto json = core::runtime_state_summary_to_json(summary, "  ");
    TEST_CHECK(json.find("\"specific_humidity\"") != std::string::npos);
    TEST_CHECK(json.find("\"cloud_water_mixing_ratio\"") != std::string::npos);
    TEST_CHECK(json.find("\"rain_water_mixing_ratio\"") != std::string::npos);
    TEST_CHECK(json.find("\"total_water_path_sum_kg_m2\"") != std::string::npos);
    TEST_CHECK(json.find("\"accumulated_surface_precipitation_sum_mm\"") !=
               std::string::npos);
    TEST_CHECK(json.find("\"mixing_ratio\"") != std::string::npos);

    return 0;
  } catch (const std::exception& ex) {
    test_fail(std::string("exception: ") + ex.what());
  }
}
