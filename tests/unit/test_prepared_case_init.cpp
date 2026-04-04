#include <cmath>
#include <vector>

#include "gwm/dycore/passive_tracer.hpp"
#include "gwm/ingest/prepared_case_init.hpp"

#include "test_assert.hpp"

namespace {

gwm::ingest::AnalysisStateIR make_analysis() {
  using namespace gwm;
  using namespace gwm::ingest;

  AnalysisStateIR analysis{};
  analysis.source = SourceKind::HRRR;
  analysis.grid.nx = 2;
  analysis.grid.ny = 2;
  analysis.grid.nz = 2;
  analysis.grid.dx = 3000.0;
  analysis.grid.dy = 3000.0;
  analysis.grid.z_top = 2000.0;
  analysis.cycle_time_utc = "2026-04-03T00:00:00Z";
  analysis.valid_time_utc = "2026-04-03T00:00:00Z";
  analysis.forecast_offset_seconds = 0;
  analysis.metadata["status"] = "populated";

  const auto count3d = analysis.grid.cell_count_3d();
  const auto count2d = analysis.grid.cell_count_2d();
  analysis.atmosphere.values["u_wind"] = std::vector<real>(count3d, 1.0f);
  analysis.atmosphere.values["v_wind"] = std::vector<real>(count3d, 2.0f);
  analysis.atmosphere.values["w_wind"] = std::vector<real>(count3d, 0.0f);
  analysis.atmosphere.values["air_temperature"] =
      std::vector<real>(count3d, 300.0f);
  analysis.atmosphere.values["specific_humidity"] =
      std::vector<real>(count3d, 0.005f);
  analysis.atmosphere.values["air_pressure"] =
      std::vector<real>(count3d, 90000.0f);
  analysis.atmosphere.values["geopotential_height"] =
      std::vector<real>(count3d, 1000.0f);
  analysis.surface.values["surface_pressure"] =
      std::vector<real>(count2d, 95000.0f);
  analysis.surface.values["air_temperature_2m"] =
      std::vector<real>(count2d, 299.0f);
  analysis.surface.values["specific_humidity_2m"] =
      std::vector<real>(count2d, 0.006f);
  analysis.surface.values["u_wind_10m"] = std::vector<real>(count2d, 5.0f);
  analysis.surface.values["v_wind_10m"] = std::vector<real>(count2d, 1.0f);
  analysis.surface.values["skin_temperature"] =
      std::vector<real>(count2d, 301.0f);
  analysis.static_surface.values["terrain_height"] = {100.0f, 120.0f, 140.0f,
                                                      160.0f};
  analysis.static_surface.values["land_mask"] =
      std::vector<real>(count2d, 1.0f);
  analysis.static_surface.values["land_use_index"] =
      std::vector<real>(count2d, 7.0f);
  return analysis;
}

gwm::ingest::BoundaryCacheIR make_boundary_cache() {
  using namespace gwm;
  using namespace gwm::ingest;

  const auto analysis = make_analysis();

  BoundaryCacheIR cache{};
  cache.source = analysis.source;
  cache.grid = analysis.grid;
  cache.cycle_time_utc = analysis.cycle_time_utc;
  cache.boundary_interval_seconds = 3600;
  cache.metadata["status"] = "populated";

  BoundarySnapshotIR t0{};
  t0.forecast_offset_seconds = 0;
  t0.valid_time_utc = analysis.valid_time_utc;
  t0.atmosphere = analysis.atmosphere;
  t0.surface = analysis.surface;

  BoundarySnapshotIR t1 = t0;
  t1.forecast_offset_seconds = 3600;
  t1.valid_time_utc = "2026-04-03T01:00:00Z";
  t1.atmosphere.values["u_wind"] =
      std::vector<real>(analysis.grid.cell_count_3d(), 3.0f);

  cache.snapshots.push_back(t0);
  cache.snapshots.push_back(t1);
  return cache;
}

}  // namespace

int main() {
  using namespace gwm;
  using namespace gwm::ingest;

  constexpr real kReferencePressure = 100000.0f;
  constexpr real kDryRd = 287.05f;
  constexpr real kKappa = 0.2854f;

  const auto analysis = make_analysis();
  PreparedCaseInitConfig config{};
  const auto layout = make_prepared_case_layout(analysis, config);
  TEST_CHECK(layout.size() == 1);
  TEST_CHECK(layout.front().nx_local() == analysis.grid.nx);
  TEST_CHECK(layout.front().ny_local() == analysis.grid.ny);

  const auto metrics = make_prepared_case_metrics(analysis, config);
  TEST_NEAR(metrics.terrain_height(0, 0), 100.0f, 1.0e-6f);
  TEST_NEAR(metrics.terrain_height(1, 1), 160.0f, 1.0e-6f);

  auto states = make_dry_states_from_analysis(analysis, layout);
  TEST_CHECK(states.size() == 1);
  const real expected_rho = 90000.0f / (kDryRd * 300.0f);
  const real expected_theta =
      300.0f * std::pow(kReferencePressure / 90000.0f, kKappa);
  TEST_NEAR(states[0].rho_d(0, 0, 0), expected_rho, 1.0e-5f);
  TEST_NEAR(states[0].rho_theta_m(0, 0, 0) / states[0].rho_d(0, 0, 0),
            expected_theta, 1.0e-3f);
  TEST_NEAR(states[0].mom_u.storage()(0, 0, 0), expected_rho, 1.0e-5f);
  TEST_NEAR(states[0].mom_v.storage()(0, 0, 0), expected_rho * 2.0f, 1.0e-5f);

  auto tracers = dycore::make_specific_humidity_tracers_from_global_field(
      states, layout, analysis.grid.nx, analysis.grid.ny,
      analysis.atmosphere.values.at("specific_humidity"), "prepared_case_qv");
  TEST_CHECK(tracers.size() == 1);
  TEST_NEAR(
      tracers[0].mass(gwm::state::kSpecificHumidityTracerName)(0, 0, 0),
      expected_rho * 0.005f, 1.0e-5f);
  TEST_NEAR(
      tracers[0].mass(gwm::state::kSpecificHumidityTracerName)(1, 1, 1),
      expected_rho * 0.005f, 1.0e-5f);

  auto warm_rain_tracers = make_warm_rain_tracers_from_analysis(
      analysis, states, layout, "prepared_case_warm_rain");
  TEST_CHECK(warm_rain_tracers.size() == 1);
  TEST_CHECK(warm_rain_tracers[0].size() == 3);
  TEST_NEAR(
      warm_rain_tracers[0].mass(gwm::state::kSpecificHumidityTracerName)(0, 0, 0),
      expected_rho * 0.005f, 1.0e-5f);
  TEST_NEAR(
      warm_rain_tracers[0].mass(gwm::state::kCloudWaterTracerName)(0, 0, 0),
      0.0f, 1.0e-6f);
  TEST_NEAR(
      warm_rain_tracers[0].mass(gwm::state::kRainWaterTracerName)(0, 0, 0),
      0.0f, 1.0e-6f);

  states[0].fill_constant(0.5f, 290.0f, 0.0f, 0.0f, 0.0f);
  PreparedCaseBoundaryUpdater boundary_updater(make_analysis(),
                                               make_boundary_cache(), config);
  boundary_updater.set_step_start_time(1800.0f);
  boundary_updater.apply(states, layout, 0.0f);
  TEST_NEAR(states[0].mom_u.storage()(0, 0, 0), expected_rho * 2.0f, 1.0e-4f);

  return 0;
}
