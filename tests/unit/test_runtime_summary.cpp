#include <cstdlib>
#include <exception>
#include <string>
#include <vector>

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

    std::vector<state::TracerState> tracers;
    tracers.emplace_back(state::make_warm_rain_registry(), 2, 1, 1, 1,
                         "runtime_summary");
    tracers.front().fill_zero();
    tracers.front().mass(state::kSpecificHumidityTracerName)(0, 0, 0) = 0.010f;
    tracers.front().mass(state::kSpecificHumidityTracerName)(1, 0, 0) = 0.012f;
    tracers.front().mass(state::kCloudWaterTracerName)(0, 0, 0) = 0.001f;
    tracers.front().mass(state::kRainWaterTracerName)(1, 0, 0) = 0.002f;

    const auto summary = core::summarize_runtime_state(states, tracers);
    const auto& qv = require_tracer(summary, state::kSpecificHumidityTracerName);
    const auto& qc = require_tracer(summary, state::kCloudWaterTracerName);
    const auto& qr = require_tracer(summary, state::kRainWaterTracerName);
    TEST_NEAR(summary.moisture.vapor_water_mass, 0.022, 1.0e-6);
    TEST_NEAR(summary.moisture.cloud_water_mass, 0.001, 1.0e-6);
    TEST_NEAR(summary.moisture.rain_water_mass, 0.002, 1.0e-6);
    TEST_NEAR(summary.moisture.condensed_water_mass, 0.003, 1.0e-6);
    TEST_NEAR(summary.moisture.total_water_mass, 0.025, 1.0e-6);
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
    TEST_CHECK(json.find("\"total_water_mass\"") != std::string::npos);
    TEST_CHECK(json.find("\"mixing_ratio\"") != std::string::npos);

    return 0;
  } catch (const std::exception& ex) {
    test_fail(std::string("exception: ") + ex.what());
  }
}
