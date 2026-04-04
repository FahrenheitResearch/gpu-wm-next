#include <cstdlib>
#include <exception>
#include <string>
#include <vector>

#include "gwm/core/runtime_summary.hpp"
#include "gwm/physics/warm_rain.hpp"
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
    states.emplace_back(3, 2, 2, 1, "warm_rain_summary");
    states.front().fill_constant(1.0f, 298.0f, 0.0f, 0.0f, 0.0f);

    std::vector<state::TracerState> tracers;
    tracers.emplace_back(state::make_warm_rain_registry(), 3, 2, 2, 1,
                         "warm_rain_summary");
    tracers.front().fill_zero();
    for (int k = 0; k < 2; ++k) {
      for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 3; ++i) {
          tracers.front().mass(state::kSpecificHumidityTracerName)(i, j, k) =
              (i == 1 ? 0.030f : 0.020f);
        }
      }
    }

    const auto before = core::summarize_runtime_state(states, tracers);

    physics::WarmRainConfig config{};
    config.dt = 2.0f;
    config.cloud_autoconversion_rate = 0.20f;
    config.rain_evaporation_rate = 0.0f;
    physics::apply_warm_rain_microphysics(states, tracers, config);
    physics::apply_warm_rain_microphysics(states, tracers, config);

    const auto after = core::summarize_runtime_state(states, tracers);
    const auto& qv = require_tracer(after, state::kSpecificHumidityTracerName);
    const auto& qc = require_tracer(after, state::kCloudWaterTracerName);
    const auto& qr = require_tracer(after, state::kRainWaterTracerName);
    TEST_NEAR(after.moisture.total_water_mass, before.moisture.total_water_mass,
              1.0e-4);
    TEST_CHECK(after.dry.total_rho_theta_m > before.dry.total_rho_theta_m);
    TEST_CHECK(after.moisture.condensed_water_mass > 0.0);
    TEST_CHECK(after.moisture.rain_water_mass > 0.0);
    TEST_CHECK(qv.mixing_ratio.min >= -1.0e-8);
    TEST_CHECK(qc.mixing_ratio.min >= -1.0e-8);
    TEST_CHECK(qr.mixing_ratio.min >= -1.0e-8);
    TEST_CHECK(after.moisture.condensed_water_mass >=
               after.moisture.rain_water_mass);
    TEST_CHECK(qr.mixing_ratio.max > 0.0);

    const auto json = core::runtime_state_summary_to_json(after, "  ");
    TEST_CHECK(json.find("\"moisture\"") != std::string::npos);
    TEST_CHECK(json.find("\"rain_water_mass\"") != std::string::npos);
    TEST_CHECK(json.find("\"tracers\"") != std::string::npos);

    return 0;
  } catch (const std::exception& ex) {
    test_fail(std::string("exception: ") + ex.what());
  }
}
