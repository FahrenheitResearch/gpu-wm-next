#include <cmath>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/domain/grid_metrics.hpp"
#include "gwm/physics/warm_rain.hpp"
#include "gwm/state/tracer_registry.hpp"

#include "test_assert.hpp"

namespace {

gwm::dycore::DryState make_single_cell_state(gwm::real rho, gwm::real theta) {
  gwm::dycore::DryState state(1, 1, 1, 1, "warm_rain_test");
  state.fill_constant(rho, theta, 0.0f, 0.0f, 0.0f);
  return state;
}

gwm::state::TracerState make_warm_rain_tracers(gwm::real rho, gwm::real qv,
                                               gwm::real qc = 0.0f,
                                               gwm::real qr = 0.0f) {
  gwm::state::TracerState tracers(gwm::state::make_warm_rain_registry(), 1, 1, 1,
                                  1, "warm_rain_test");
  tracers.fill_zero();
  tracers.mass(gwm::state::kSpecificHumidityTracerName)(0, 0, 0) = rho * qv;
  tracers.mass(gwm::state::kCloudWaterTracerName)(0, 0, 0) = rho * qc;
  tracers.mass(gwm::state::kRainWaterTracerName)(0, 0, 0) = rho * qr;
  return tracers;
}

gwm::real mixing_ratio(const gwm::state::TracerState& tracers,
                       const char* name) {
  return tracers.mass(name)(0, 0, 0);
}

}  // namespace

int main() {
  try {
    using namespace gwm;

    {
      std::vector<dycore::DryState> states;
      states.push_back(make_single_cell_state(1.0f, 300.0f));
      std::vector<state::TracerState> tracers;
      tracers.push_back(make_warm_rain_tracers(1.0f, 0.030f, 0.0f, 0.0f));

      const real theta_before =
          states.front().rho_theta_m(0, 0, 0) / states.front().rho_d(0, 0, 0);
      const real water_before =
          mixing_ratio(tracers.front(), state::kSpecificHumidityTracerName) +
          mixing_ratio(tracers.front(), state::kCloudWaterTracerName) +
          mixing_ratio(tracers.front(), state::kRainWaterTracerName);

      physics::WarmRainConfig config{};
      config.dt = 2.0f;
      config.cloud_autoconversion_rate = 0.0f;
      physics::apply_warm_rain_microphysics(states, tracers, config);

      const real qv =
          mixing_ratio(tracers.front(), state::kSpecificHumidityTracerName);
      const real qc = mixing_ratio(tracers.front(), state::kCloudWaterTracerName);
      const real qr = mixing_ratio(tracers.front(), state::kRainWaterTracerName);
      const real theta_after =
          states.front().rho_theta_m(0, 0, 0) / states.front().rho_d(0, 0, 0);
      const real water_after = qv + qc + qr;

      TEST_CHECK(qv < 0.030f);
      TEST_CHECK(qc > 0.0f);
      TEST_NEAR(qr, 0.0f, 1.0e-6f);
      TEST_NEAR(water_after, water_before, 1.0e-5f);
      TEST_CHECK(theta_after > theta_before);
    }

    {
      std::vector<dycore::DryState> states;
      states.emplace_back(1, 1, 2, 1, "warm_rain_nominal_fallout");
      states.front().fill_constant(1.0f, 295.0f, 0.0f, 0.0f, 0.0f);

      std::vector<state::TracerState> tracers;
      tracers.emplace_back(state::make_warm_rain_registry(), 1, 1, 2, 1,
                           "warm_rain_nominal_fallout");
      tracers.front().fill_zero();
      tracers.front().mass(state::kRainWaterTracerName)(0, 0, 0) = 0.004f;
      tracers.front().mass(state::kRainWaterTracerName)(0, 0, 1) = 0.010f;

      std::vector<physics::WarmRainSurfaceAccumulation> accumulations(1);
      accumulations.front().reset(1, 1);

      physics::WarmRainConfig config{};
      config.dt = 1.0f;
      config.cloud_autoconversion_rate = 0.0f;
      config.rain_evaporation_rate = 0.0f;
      config.rain_terminal_velocity = 250.0f;
      config.sedimentation_layer_depth_m = 1000.0f;
      config.enable_latent_heating = false;
      physics::apply_warm_rain_microphysics(states, tracers, config,
                                            &accumulations);

      TEST_NEAR(tracers.front().mass(state::kRainWaterTracerName)(0, 0, 1),
                0.0075f, 1.0e-6f);
      TEST_NEAR(tracers.front().mass(state::kRainWaterTracerName)(0, 0, 0),
                0.0055f, 1.0e-6f);
      TEST_NEAR(accumulations.front().liquid_precipitation_kg_m2(0, 0, 0), 1.0f,
                1.0e-6f);
    }

    {
      std::vector<dycore::DryState> states;
      states.push_back(make_single_cell_state(1.0f, 295.0f));
      std::vector<state::TracerState> tracers;
      tracers.push_back(make_warm_rain_tracers(1.0f, 0.006f, 0.002f, 0.001f));

      const real water_before =
          mixing_ratio(tracers.front(), state::kSpecificHumidityTracerName) +
          mixing_ratio(tracers.front(), state::kCloudWaterTracerName) +
          mixing_ratio(tracers.front(), state::kRainWaterTracerName);

      physics::WarmRainConfig config{};
      config.dt = 2.0f;
      config.cloud_autoconversion_threshold = 1.0e-4f;
      config.cloud_autoconversion_rate = 1.0f;
      config.rain_evaporation_rate = 0.0f;
      config.enable_latent_heating = false;
      physics::apply_warm_rain_microphysics(states, tracers, config);

      const real qc = mixing_ratio(tracers.front(), state::kCloudWaterTracerName);
      const real qr = mixing_ratio(tracers.front(), state::kRainWaterTracerName);
      const real water_after =
          mixing_ratio(tracers.front(), state::kSpecificHumidityTracerName) + qc + qr;

      TEST_CHECK(qc < 0.002f);
      TEST_CHECK(qr > 0.001f);
      TEST_NEAR(water_after, water_before, 1.0e-5f);
    }

    {
      std::vector<dycore::DryState> states;
      states.push_back(make_single_cell_state(1.0f, 295.0f));
      std::vector<state::TracerState> tracers;
      tracers.push_back(make_warm_rain_tracers(1.0f, 0.0f, 0.0f, 0.010f));
      std::vector<physics::WarmRainSurfaceAccumulation> accumulations(1);
      accumulations.front().reset(1, 1);
      const auto layout =
          comm::VirtualRankLayout::build(1, 1, 1, 1, 1, 1);
      const auto metrics = domain::GridMetrics::make_hybrid_height(
          1, 1, 1, 1000.0, 1000.0, 1000.0, std::vector<real>{0.0f});

      physics::WarmRainConfig config{};
      config.dt = 1.0f;
      config.cloud_autoconversion_rate = 0.0f;
      config.rain_evaporation_rate = 0.0f;
      config.rain_terminal_velocity = 1500.0f;
      config.enable_latent_heating = false;
      physics::apply_warm_rain_microphysics(states, tracers, layout, metrics,
                                            config, &accumulations);

      TEST_NEAR(mixing_ratio(tracers.front(), state::kRainWaterTracerName), 0.0f,
                1.0e-6f);
      TEST_NEAR(accumulations.front().liquid_precipitation_kg_m2(0, 0, 0), 10.0f,
                1.0e-6f);
    }

    return 0;
  } catch (const std::exception& ex) {
    test_fail(std::string("exception: ") + ex.what());
  }
}
