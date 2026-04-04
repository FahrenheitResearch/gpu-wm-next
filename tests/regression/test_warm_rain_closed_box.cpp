#include <vector>

#include "gwm/physics/warm_rain.hpp"
#include "gwm/state/tracer_registry.hpp"

#include "test_assert.hpp"

int main() {
  try {
    using namespace gwm;

    std::vector<dycore::DryState> states;
    states.emplace_back(4, 3, 2, 1, "warm_rain_box");
    states.front().fill_constant(1.0f, 300.0f, 0.0f, 0.0f, 0.0f);

    std::vector<state::TracerState> tracers;
    tracers.emplace_back(state::make_warm_rain_registry(), 4, 3, 2, 1,
                         "warm_rain_box");
    tracers.front().fill_zero();
    for (int k = 0; k < 2; ++k) {
      for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 4; ++i) {
          tracers.front().mass(state::kSpecificHumidityTracerName)(i, j, k) =
              (j == 1 ? 0.030f : 0.020f);
          tracers.front().mass(state::kCloudWaterTracerName)(i, j, k) = 0.0f;
          tracers.front().mass(state::kRainWaterTracerName)(i, j, k) = 0.0f;
        }
      }
    }

    const double water_before =
        tracers.front().total_mass(
            tracers.front().registry().find(state::kSpecificHumidityTracerName)
                .value()) +
        tracers.front().total_mass(
            tracers.front().registry().find(state::kCloudWaterTracerName)
                .value()) +
        tracers.front().total_mass(
            tracers.front().registry().find(state::kRainWaterTracerName).value());
    const double theta_before = states.front().total_rho_theta_m();

    physics::WarmRainConfig config{};
    config.dt = 2.0f;
    config.cloud_autoconversion_rate = 0.20f;
    config.rain_evaporation_rate = 0.0f;
    physics::apply_warm_rain_microphysics(states, tracers, config);
    physics::apply_warm_rain_microphysics(states, tracers, config);

    const double water_after =
        tracers.front().total_mass(
            tracers.front().registry().find(state::kSpecificHumidityTracerName)
                .value()) +
        tracers.front().total_mass(
            tracers.front().registry().find(state::kCloudWaterTracerName)
                .value()) +
        tracers.front().total_mass(
            tracers.front().registry().find(state::kRainWaterTracerName).value());
    const double theta_after = states.front().total_rho_theta_m();

    TEST_NEAR(water_after, water_before, 1.0e-4);
    TEST_CHECK(theta_after > theta_before);
    const double condensed_after =
        tracers.front().total_mass(
            tracers.front().registry().find(state::kCloudWaterTracerName)
                .value()) +
        tracers.front().total_mass(
            tracers.front().registry().find(state::kRainWaterTracerName)
                .value());
    TEST_CHECK(condensed_after > 0.0);
    TEST_CHECK(tracers.front().total_mass(
                   tracers.front().registry().find(state::kRainWaterTracerName)
                       .value()) > 0.0);

    return 0;
  } catch (const std::exception& ex) {
    test_fail(std::string("exception: ") + ex.what());
  }
}
