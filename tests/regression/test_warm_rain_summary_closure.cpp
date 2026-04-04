#include <cstdlib>
#include <exception>
#include <string>
#include <vector>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/core/runtime_summary.hpp"
#include "gwm/domain/grid_metrics.hpp"
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
    std::vector<physics::WarmRainSurfaceAccumulation> accumulations(1);
    accumulations.front().reset(3, 2);
    const auto layout =
        comm::VirtualRankLayout::build(3, 2, 2, 1, 1, 1);
    const auto metrics = domain::GridMetrics::make_hybrid_height(
        3, 2, 2, 1000.0, 1000.0, 2000.0,
        std::vector<real>(static_cast<std::size_t>(3 * 2), 0.0f));
    const auto before =
        core::summarize_runtime_state(states, tracers, &accumulations);
    const auto total_before =
        total_column_water_per_area(tracers.front(), layout.front(), metrics);

    physics::WarmRainConfig config{};
    config.dt = 2.0f;
    config.cloud_autoconversion_rate = 0.20f;
    config.rain_evaporation_rate = 0.0f;
    config.rain_terminal_velocity = 2.0f;
    physics::apply_warm_rain_microphysics(states, tracers, layout, metrics,
                                          config, &accumulations);
    physics::apply_warm_rain_microphysics(states, tracers, layout, metrics,
                                          config, &accumulations);

    const auto after =
        core::summarize_runtime_state(states, tracers, &accumulations);
    const auto& qv = require_tracer(after, state::kSpecificHumidityTracerName);
    const auto& qc = require_tracer(after, state::kCloudWaterTracerName);
    const auto& qr = require_tracer(after, state::kRainWaterTracerName);
    const auto total_after =
        total_column_water_per_area(tracers.front(), layout.front(), metrics) +
        accumulations.front().liquid_precipitation_kg_m2.owned_sum();
    TEST_NEAR(total_after, total_before, 1.0e-3);
    TEST_CHECK(after.dry.total_rho_theta_m > before.dry.total_rho_theta_m);
    TEST_CHECK(after.moisture.condensed_water_mass > 0.0);
    TEST_CHECK(after.moisture.accumulated_surface_precipitation_sum_mm > 0.0);
    TEST_CHECK(qv.mixing_ratio.min >= -1.0e-8);
    TEST_CHECK(qc.mixing_ratio.min >= -1.0e-8);
    TEST_CHECK(qr.mixing_ratio.min >= -1.0e-8);
    TEST_CHECK(after.moisture.max_surface_precipitation_mm > 0.0);
    TEST_CHECK(after.moisture.mean_surface_precipitation_mm > 0.0);

    const auto json = core::runtime_state_summary_to_json(after, "  ");
    TEST_CHECK(json.find("\"moisture\"") != std::string::npos);
    TEST_CHECK(json.find("\"accumulated_surface_precipitation_sum_mm\"") !=
               std::string::npos);
    TEST_CHECK(json.find("\"tracers\"") != std::string::npos);

    return 0;
  } catch (const std::exception& ex) {
    test_fail(std::string("exception: ") + ex.what());
  }
}
