#include <vector>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/domain/grid_metrics.hpp"
#include "gwm/physics/warm_rain.hpp"
#include "gwm/state/tracer_registry.hpp"

#include "test_assert.hpp"

namespace {

double total_column_water_per_area(
    const gwm::state::TracerState& tracers,
    const gwm::domain::SubdomainDescriptor& desc,
    const gwm::domain::GridMetrics& metrics) {
  const auto qv_index =
      tracers.registry().find(gwm::state::kSpecificHumidityTracerName).value();
  const auto qc_index =
      tracers.registry().find(gwm::state::kCloudWaterTracerName).value();
  const auto qr_index =
      tracers.registry().find(gwm::state::kRainWaterTracerName).value();

  double total = 0.0;
  for (int k = 0; k < desc.nz; ++k) {
    for (int j = 0; j < desc.ny_local(); ++j) {
      for (int i = 0; i < desc.nx_local(); ++i) {
        const double dz = 1.0 / static_cast<double>(
            metrics.inv_dz_cell(desc.i_begin + i, desc.j_begin + j, k));
        total += (static_cast<double>(tracers.mass(qv_index)(i, j, k)) +
                  static_cast<double>(tracers.mass(qc_index)(i, j, k)) +
                  static_cast<double>(tracers.mass(qr_index)(i, j, k))) *
                 dz;
      }
    }
  }
  return total;
}

}  // namespace

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

    std::vector<dycore::DryState> fallout_states;
    fallout_states.emplace_back(2, 1, 2, 1, "warm_rain_fallout_box");
    fallout_states.front().fill_constant(1.0f, 300.0f, 0.0f, 0.0f, 0.0f);

    std::vector<state::TracerState> fallout_tracers;
    fallout_tracers.emplace_back(state::make_warm_rain_registry(), 2, 1, 2, 1,
                                 "warm_rain_fallout_box");
    fallout_tracers.front().fill_zero();
    for (int k = 0; k < 2; ++k) {
      for (int i = 0; i < 2; ++i) {
        fallout_tracers.front().mass(state::kRainWaterTracerName)(i, 0, k) =
            (k == 1 ? 0.010f : 0.004f);
      }
    }
    std::vector<physics::WarmRainSurfaceAccumulation> accumulations(1);
    accumulations.front().reset(2, 1);
    const auto layout =
        comm::VirtualRankLayout::build(2, 1, 2, 1, 1, 1);
    const auto metrics = domain::GridMetrics::make_hybrid_height(
        2, 1, 2, 1000.0, 1000.0, 2000.0, std::vector<real>{0.0f, 0.0f});

    const auto total_before =
        total_column_water_per_area(fallout_tracers.front(), layout.front(),
                                    metrics);

    physics::WarmRainConfig fallout_config{};
    fallout_config.dt = 1.0f;
    fallout_config.condensation_relaxation = 0.0f;
    fallout_config.cloud_autoconversion_rate = 0.0f;
    fallout_config.rain_evaporation_rate = 0.0f;
    fallout_config.rain_terminal_velocity = 1000.0f;
    fallout_config.enable_latent_heating = false;
    physics::apply_warm_rain_microphysics(fallout_states, fallout_tracers,
                                          layout, metrics, fallout_config,
                                          &accumulations);

    const auto total_after_atmosphere =
        total_column_water_per_area(fallout_tracers.front(), layout.front(),
                                    metrics);
    const auto total_surface =
        accumulations.front().liquid_precipitation_kg_m2.owned_sum();
    TEST_NEAR(total_after_atmosphere + total_surface, total_before, 1.0e-5);
    TEST_CHECK(total_surface > 0.0);

    return 0;
  } catch (const std::exception& ex) {
    test_fail(std::string("exception: ") + ex.what());
  }
}
