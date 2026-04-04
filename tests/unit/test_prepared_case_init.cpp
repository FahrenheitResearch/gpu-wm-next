#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "gwm/core/dry_thermo.hpp"
#include "gwm/ingest/prepared_case_init.hpp"

#include "test_assert.hpp"

namespace {

template <typename Fn>
void fill_3d(std::vector<gwm::real>& field, int nx, int ny, int nz, Fn&& fn) {
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        field[(static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
               static_cast<std::size_t>(j)) *
                  static_cast<std::size_t>(nx) +
              static_cast<std::size_t>(i)] = fn(i, j, k);
      }
    }
  }
}

std::size_t linear_index_3d(int i, int j, int k, int nx, int ny) {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

gwm::real source_theta(const gwm::ingest::AnalysisStateIR& analysis, int i,
                       int j, int k) {
  const auto& pressure = analysis.atmosphere.values.at("air_pressure");
  const auto& temperature = analysis.atmosphere.values.at("air_temperature");
  const auto idx = linear_index_3d(i, j, k, analysis.grid.nx, analysis.grid.ny);
  return temperature[idx] *
         std::pow(gwm::core::kReferencePressure /
                      std::max(pressure[idx], 1.0f),
                  gwm::core::kKappa);
}

struct ProjectedColumn {
  std::vector<gwm::real> pressure;
  std::vector<gwm::real> rho;
};

ProjectedColumn project_column(const gwm::ingest::AnalysisStateIR& analysis,
                               const gwm::domain::GridMetrics& metrics, int i,
                               int j) {
  constexpr gwm::real kGravity = 9.81f;

  ProjectedColumn column{};
  column.pressure.resize(static_cast<std::size_t>(analysis.grid.nz));
  column.rho.resize(static_cast<std::size_t>(analysis.grid.nz));

  const auto idx0 = linear_index_3d(i, j, 0, analysis.grid.nx, analysis.grid.ny);
  column.pressure[0] =
      std::max(analysis.atmosphere.values.at("air_pressure")[idx0], 1.0f);
  column.rho[0] = gwm::core::dry_rho_from_pressure_theta(
      column.pressure[0], source_theta(analysis, i, j, 0));

  for (int k = 1; k < analysis.grid.nz; ++k) {
    const gwm::real theta_k = source_theta(analysis, i, j, k);
    const gwm::real dz =
        metrics.z_center(i, j, k) - metrics.z_center(i, j, k - 1);
    gwm::real p_next = std::max<gwm::real>(
        1.0f, column.pressure[static_cast<std::size_t>(k - 1)] -
                  kGravity * column.rho[static_cast<std::size_t>(k - 1)] * dz);
    for (int iter = 0; iter < 4; ++iter) {
      const gwm::real rho_next =
          gwm::core::dry_rho_from_pressure_theta(p_next, theta_k);
      p_next = column.pressure[static_cast<std::size_t>(k - 1)] -
               kGravity * 0.5f *
                   (column.rho[static_cast<std::size_t>(k - 1)] + rho_next) * dz;
      p_next = std::max(p_next, 1.0f);
    }
    column.pressure[static_cast<std::size_t>(k)] = p_next;
    column.rho[static_cast<std::size_t>(k)] =
        gwm::core::dry_rho_from_pressure_theta(p_next, theta_k);
  }

  return column;
}

gwm::real expected_face_mom_u(const gwm::ingest::AnalysisStateIR& analysis,
                              const gwm::domain::GridMetrics& metrics,
                              int i_face, int j, int k) {
  const int left_i = std::clamp(i_face - 1, 0, analysis.grid.nx - 1);
  const int right_i = std::clamp(i_face, 0, analysis.grid.nx - 1);
  const auto left = project_column(analysis, metrics, left_i, j);
  const auto right = project_column(analysis, metrics, right_i, j);
  const auto idx_left =
      linear_index_3d(left_i, j, k, analysis.grid.nx, analysis.grid.ny);
  const auto idx_right =
      linear_index_3d(right_i, j, k, analysis.grid.nx, analysis.grid.ny);
  const gwm::real rho_face =
      0.5f * (left.rho[static_cast<std::size_t>(k)] +
              right.rho[static_cast<std::size_t>(k)]);
  const gwm::real u_face =
      0.5f * (analysis.atmosphere.values.at("u_wind")[idx_left] +
              analysis.atmosphere.values.at("u_wind")[idx_right]);
  return rho_face * u_face;
}

gwm::real expected_face_mom_v(const gwm::ingest::AnalysisStateIR& analysis,
                              const gwm::domain::GridMetrics& metrics,
                              int i, int j_face, int k) {
  const int south_j = std::clamp(j_face - 1, 0, analysis.grid.ny - 1);
  const int north_j = std::clamp(j_face, 0, analysis.grid.ny - 1);
  const auto south = project_column(analysis, metrics, i, south_j);
  const auto north = project_column(analysis, metrics, i, north_j);
  const auto idx_south =
      linear_index_3d(i, south_j, k, analysis.grid.nx, analysis.grid.ny);
  const auto idx_north =
      linear_index_3d(i, north_j, k, analysis.grid.nx, analysis.grid.ny);
  const gwm::real rho_face =
      0.5f * (south.rho[static_cast<std::size_t>(k)] +
              north.rho[static_cast<std::size_t>(k)]);
  const gwm::real v_face =
      0.5f * (analysis.atmosphere.values.at("v_wind")[idx_south] +
              analysis.atmosphere.values.at("v_wind")[idx_north]);
  return rho_face * v_face;
}

std::vector<gwm::real> expected_projected_mom_w_faces(
    const gwm::ingest::AnalysisStateIR& analysis,
    const gwm::domain::GridMetrics& metrics, int i, int j) {
  std::vector<gwm::real> faces(static_cast<std::size_t>(analysis.grid.nz + 1),
                               0.0f);
  gwm::real column_mass_bias = 0.0f;
  gwm::real column_height = 0.0f;

  for (int k = 0; k < analysis.grid.nz; ++k) {
    const gwm::real hdiv =
        (expected_face_mom_u(analysis, metrics, i + 1, j, k) -
         expected_face_mom_u(analysis, metrics, i, j, k)) /
            static_cast<gwm::real>(analysis.grid.dx) +
        (expected_face_mom_v(analysis, metrics, i, j + 1, k) -
         expected_face_mom_v(analysis, metrics, i, j, k)) /
            static_cast<gwm::real>(analysis.grid.dy);
    const gwm::real inv_dz = metrics.inv_dz_cell(i, j, k);
    const gwm::real dz = 1.0f / std::max(inv_dz, 1.0e-6f);
    column_mass_bias += hdiv * dz;
    column_height += dz;
  }

  const gwm::real mean_divergence =
      column_mass_bias / std::max(column_height, 1.0e-6f);
  for (int k = 0; k < analysis.grid.nz; ++k) {
    const gwm::real hdiv =
        (expected_face_mom_u(analysis, metrics, i + 1, j, k) -
         expected_face_mom_u(analysis, metrics, i, j, k)) /
            static_cast<gwm::real>(analysis.grid.dx) +
        (expected_face_mom_v(analysis, metrics, i, j + 1, k) -
         expected_face_mom_v(analysis, metrics, i, j, k)) /
            static_cast<gwm::real>(analysis.grid.dy);
    const gwm::real inv_dz = metrics.inv_dz_cell(i, j, k);
    const gwm::real dz = 1.0f / std::max(inv_dz, 1.0e-6f);
    faces[static_cast<std::size_t>(k + 1)] =
        faces[static_cast<std::size_t>(k)] - (hdiv - mean_divergence) * dz;
  }
  faces.back() = 0.0f;
  return faces;
}

void apply_constant_offset(std::vector<gwm::real>& values, gwm::real delta) {
  for (auto& value : values) {
    value += delta;
  }
}

gwm::ingest::AnalysisStateIR make_analysis() {
  using namespace gwm;
  using namespace gwm::ingest;

  AnalysisStateIR analysis{};
  analysis.source = SourceKind::HRRR;
  analysis.grid.nx = 2;
  analysis.grid.ny = 2;
  analysis.grid.nz = 3;
  analysis.grid.dx = 3000.0;
  analysis.grid.dy = 3000.0;
  analysis.grid.z_top = 3000.0;
  analysis.cycle_time_utc = "2026-04-03T00:00:00Z";
  analysis.valid_time_utc = "2026-04-03T00:00:00Z";
  analysis.forecast_offset_seconds = 0;
  analysis.metadata["status"] = "populated";

  const auto count3d = analysis.grid.cell_count_3d();
  const auto count2d = analysis.grid.cell_count_2d();
  analysis.atmosphere.values["air_pressure"] =
      std::vector<real>(count3d, 90000.0f);
  analysis.atmosphere.values["air_temperature"] =
      std::vector<real>(count3d, 300.0f);
  analysis.atmosphere.values["specific_humidity"] =
      std::vector<real>(count3d, 0.005f);
  analysis.atmosphere.values["u_wind"] = std::vector<real>(count3d, 6.0f);
  analysis.atmosphere.values["v_wind"] = std::vector<real>(count3d, -3.0f);
  analysis.atmosphere.values["w_wind"] = std::vector<real>(count3d, 8.0f);
  analysis.atmosphere.values["geopotential_height"] =
      std::vector<real>(count3d, 0.0f);
  fill_3d(analysis.atmosphere.values["geopotential_height"], analysis.grid.nx,
          analysis.grid.ny, analysis.grid.nz,
          [](int i, int j, int k) {
            return 120.0f + 700.0f * static_cast<gwm::real>(k) +
                   40.0f * static_cast<gwm::real>(i) +
                   25.0f * static_cast<gwm::real>(j);
          });

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
  analysis.static_surface.values["terrain_height"] = {0.0f, 180.0f, 320.0f,
                                                      560.0f};
  analysis.static_surface.values["land_mask"] = std::vector<real>(count2d, 1.0f);
  analysis.static_surface.values["land_use_index"] =
      std::vector<real>(count2d, 7.0f);
  return analysis;
}

gwm::ingest::AnalysisStateIR make_boundary_analysis() {
  auto analysis = make_analysis();
  apply_constant_offset(analysis.atmosphere.values["air_pressure"], -900.0f);
  apply_constant_offset(analysis.atmosphere.values["air_temperature"], 1.0f);
  apply_constant_offset(analysis.atmosphere.values["specific_humidity"],
                        0.00125f);
  apply_constant_offset(analysis.atmosphere.values["u_wind"], 1.5f);
  apply_constant_offset(analysis.atmosphere.values["v_wind"], -0.75f);
  apply_constant_offset(analysis.atmosphere.values["w_wind"], 10.0f);
  apply_constant_offset(analysis.atmosphere.values["geopotential_height"],
                        15.0f);
  return analysis;
}

gwm::ingest::AnalysisStateIR make_late_boundary_analysis() {
  auto analysis = make_analysis();
  apply_constant_offset(analysis.atmosphere.values["air_pressure"], -1800.0f);
  apply_constant_offset(analysis.atmosphere.values["air_temperature"], 2.0f);
  apply_constant_offset(analysis.atmosphere.values["specific_humidity"],
                        0.0025f);
  apply_constant_offset(analysis.atmosphere.values["u_wind"], 3.0f);
  apply_constant_offset(analysis.atmosphere.values["v_wind"], -1.5f);
  apply_constant_offset(analysis.atmosphere.values["w_wind"], 20.0f);
  apply_constant_offset(analysis.atmosphere.values["geopotential_height"],
                        30.0f);
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
  apply_constant_offset(t1.atmosphere.values["air_pressure"], -1800.0f);
  apply_constant_offset(t1.atmosphere.values["air_temperature"], 2.0f);
  apply_constant_offset(t1.atmosphere.values["specific_humidity"], 0.0025f);
  apply_constant_offset(t1.atmosphere.values["u_wind"], 3.0f);
  apply_constant_offset(t1.atmosphere.values["v_wind"], -1.5f);
  apply_constant_offset(t1.atmosphere.values["w_wind"], 20.0f);
  apply_constant_offset(t1.atmosphere.values["geopotential_height"], 30.0f);

  cache.snapshots.push_back(t0);
  cache.snapshots.push_back(t1);
  return cache;
}

void check_balance_diag(const gwm::ingest::PreparedCaseBalanceDiagnostics& diag) {
  std::cerr << "balance diag:"
            << " eos=" << diag.max_rel_eos
            << " hydro=" << diag.max_rel_hydrostatic
            << " fast_w=" << diag.max_rel_fast_vertical
            << " mass_div_abs=" << diag.max_abs_mass_divergence
            << " mass_div_rel=" << diag.max_rel_mass_divergence
            << " mom_w_bottom=" << diag.max_abs_mom_w_bottom
            << " mom_w_top=" << diag.max_abs_mom_w_top
            << " tracer_closure=" << diag.max_abs_tracer_closure;
  if (diag.max_abs_z_src_minus_metric.has_value()) {
    std::cerr << " z_src_metric=" << *diag.max_abs_z_src_minus_metric;
  }
  std::cerr << std::endl;
  TEST_CHECK(diag.max_rel_eos <= 1.0e-6f);
  TEST_CHECK(diag.max_rel_hydrostatic <= 1.0e-6f);
  TEST_CHECK(diag.max_rel_fast_vertical <= 1.0e-6f);
  TEST_CHECK(diag.max_abs_mass_divergence <= 1.0e-6f);
  TEST_CHECK(diag.max_rel_mass_divergence <= 1.0e-3f);
  TEST_CHECK(diag.max_abs_mom_w_bottom == 0.0f);
  TEST_CHECK(diag.max_abs_mom_w_top == 0.0f);
  TEST_CHECK(diag.max_abs_tracer_closure <= 1.0e-6f);
  TEST_CHECK(diag.max_abs_z_src_minus_metric.has_value());
}

}  // namespace

int main() {
  using namespace gwm;
  using namespace gwm::ingest;

  const auto analysis = make_analysis();
  PreparedCaseInitConfig config{};
  const auto layout = make_prepared_case_layout(analysis, config);
  TEST_CHECK(layout.size() == 1);
  TEST_CHECK(layout.front().nx_local() == analysis.grid.nx);
  TEST_CHECK(layout.front().ny_local() == analysis.grid.ny);

  const auto metrics = make_prepared_case_metrics(analysis, config);
  TEST_NEAR(metrics.terrain_height(0, 0), 0.0f, 1.0e-6f);
  TEST_NEAR(metrics.terrain_height(1, 1), 560.0f, 1.0e-6f);

  auto states = make_dry_states_from_analysis(analysis, metrics, layout);
  TEST_CHECK(states.size() == 1);

  auto warm_rain_tracers = make_warm_rain_tracers_from_analysis(
      analysis, states, layout, "prepared_case_warm_rain");
  TEST_CHECK(warm_rain_tracers.size() == 1);
  TEST_CHECK(warm_rain_tracers[0].size() == 3);

  const auto startup_diag = diagnose_prepared_case_balance(
      analysis, metrics, states, layout, &warm_rain_tracers);
  check_balance_diag(startup_diag);

  const auto startup_column_00 = project_column(analysis, metrics, 0, 0);
  const auto startup_column_11 = project_column(analysis, metrics, 1, 1);
  TEST_NEAR(states[0].rho_d(0, 0, 0), startup_column_00.rho[0], 1.0e-5f);
  TEST_NEAR(states[0].rho_d(1, 1, 2), startup_column_11.rho[2], 1.0e-5f);
  TEST_NEAR(states[0].rho_theta_m(0, 0, 0) / states[0].rho_d(0, 0, 0),
            source_theta(analysis, 0, 0, 0), 1.0e-5f);
  TEST_NEAR(states[0].rho_theta_m(1, 1, 2) / states[0].rho_d(1, 1, 2),
            source_theta(analysis, 1, 1, 2), 1.0e-5f);
  auto analysis_with_shifted_w = analysis;
  fill_3d(analysis_with_shifted_w.atmosphere.values["w_wind"], analysis.grid.nx,
          analysis.grid.ny, analysis.grid.nz,
          [](int i, int j, int k) {
            return 50.0f + 2.0f * static_cast<gwm::real>(i) -
                   3.0f * static_cast<gwm::real>(j) +
                   5.0f * static_cast<gwm::real>(k);
          });
  const auto states_with_shifted_w =
      make_dry_states_from_analysis(analysis_with_shifted_w, metrics, layout);
  for (int k = 0; k < analysis.grid.nz; ++k) {
    TEST_NEAR(states_with_shifted_w[0].mom_u.storage()(0, 0, k),
              states[0].mom_u.storage()(0, 0, k), 1.0e-6f);
    TEST_NEAR(states_with_shifted_w[0].mom_v.storage()(0, 0, k),
              states[0].mom_v.storage()(0, 0, k), 1.0e-6f);
  }
  for (int k_face = 0; k_face <= analysis.grid.nz; ++k_face) {
    TEST_NEAR(states_with_shifted_w[0].mom_w.storage()(0, 0, k_face),
              states[0].mom_w.storage()(0, 0, k_face), 1.0e-6f);
    TEST_NEAR(states_with_shifted_w[0].mom_w.storage()(1, 1, k_face),
              states[0].mom_w.storage()(1, 1, k_face), 1.0e-6f);
  }
  TEST_NEAR(warm_rain_tracers[0]
                .mass(gwm::state::kSpecificHumidityTracerName)(0, 0, 0),
            states[0].rho_d(0, 0, 0) *
                analysis.atmosphere.values.at("specific_humidity")
                    [linear_index_3d(0, 0, 0, analysis.grid.nx, analysis.grid.ny)],
            1.0e-5f);
  TEST_NEAR(warm_rain_tracers[0]
                .mass(gwm::state::kSpecificHumidityTracerName)(1, 1, 2),
            states[0].rho_d(1, 1, 2) *
                analysis.atmosphere.values.at("specific_humidity")
                    [linear_index_3d(1, 1, 2, analysis.grid.nx, analysis.grid.ny)],
            1.0e-5f);
  TEST_CHECK(warm_rain_tracers[0]
                 .mass(gwm::state::kCloudWaterTracerName)
                 .owned_sum() == 0.0);
  TEST_CHECK(warm_rain_tracers[0]
                 .mass(gwm::state::kRainWaterTracerName)
                 .owned_sum() == 0.0);

  auto boundary_states =
      make_dry_states_from_analysis(analysis, metrics, layout, "boundary_seed");
  boundary_states[0].fill_constant(0.5f, 290.0f, 0.0f, 0.0f, 0.0f);

  PreparedCaseBoundaryUpdater boundary_updater(make_analysis(),
                                               make_boundary_cache(), config);
  boundary_updater.set_step_start_time(1800.0f);
  boundary_updater.apply(boundary_states, layout, 0.0f);

  const auto boundary_analysis = make_boundary_analysis();
  const auto expected_boundary_states =
      make_dry_states_from_analysis(boundary_analysis, metrics, layout);
  TEST_NEAR(boundary_states[0].rho_d(0, 0, 1),
            expected_boundary_states[0].rho_d(0, 0, 1), 1.0e-5f);
  TEST_NEAR(boundary_states[0].mom_u.storage()(0, 0, 1),
            expected_boundary_states[0].mom_u.storage()(0, 0, 1),
            1.0e-5f);
  TEST_NEAR(boundary_states[0].mom_v.storage()(1, 1, 2),
            expected_boundary_states[0].mom_v.storage()(1, 1, 2),
            1.0e-5f);
  for (int k_face = 0; k_face <= boundary_analysis.grid.nz; ++k_face) {
    TEST_NEAR(boundary_states[0].mom_w.storage()(0, 0, k_face),
              expected_boundary_states[0].mom_w.storage()(0, 0, k_face),
              1.0e-5f);
  }

  boundary_states[0].fill_constant(0.25f, 285.0f, -1.0f, 2.0f, 3.0f);
  boundary_updater.apply(boundary_states, layout, 1800.0f);
  const auto late_boundary_analysis = make_late_boundary_analysis();
  const auto expected_late_boundary_states =
      make_dry_states_from_analysis(late_boundary_analysis, metrics, layout);
  TEST_NEAR(boundary_states[0].rho_d(0, 0, 1),
            expected_late_boundary_states[0].rho_d(0, 0, 1), 1.0e-5f);
  TEST_NEAR(boundary_states[0].mom_u.storage()(0, 0, 1),
            expected_late_boundary_states[0].mom_u.storage()(0, 0, 1),
            1.0e-5f);
  TEST_NEAR(boundary_states[0].mom_v.storage()(1, 1, 2),
            expected_late_boundary_states[0].mom_v.storage()(1, 1, 2),
            1.0e-5f);

  auto boundary_warm_rain_tracers = make_warm_rain_tracers_from_analysis(
      analysis, states, layout, "prepared_case_boundary_warm_rain");
  auto& qv =
      boundary_warm_rain_tracers[0].mass(gwm::state::kSpecificHumidityTracerName);
  auto& qc =
      boundary_warm_rain_tracers[0].mass(gwm::state::kCloudWaterTracerName);
  auto& qr =
      boundary_warm_rain_tracers[0].mass(gwm::state::kRainWaterTracerName);
  qv(0, 0, 0) = 0.0f;
  qc(0, 0, 0) = 0.020f;
  qr(0, 0, 0) = 0.010f;

  PreparedCaseTracerBoundaryUpdater tracer_boundary_updater(
      make_analysis(), make_boundary_cache(), config);
  tracer_boundary_updater.set_step_start_time(1800.0f);
  tracer_boundary_updater.apply(boundary_warm_rain_tracers, layout, 0.0f);

  const auto expected_boundary_tracers = make_warm_rain_tracers_from_analysis(
      boundary_analysis, expected_boundary_states, layout,
      "expected_boundary_warm_rain");
  TEST_NEAR(boundary_warm_rain_tracers[0]
                .mass(gwm::state::kSpecificHumidityTracerName)(0, 0, 0),
            expected_boundary_tracers[0]
                .mass(gwm::state::kSpecificHumidityTracerName)(0, 0, 0),
            1.0e-5f);
  TEST_NEAR(boundary_warm_rain_tracers[0]
                .mass(gwm::state::kSpecificHumidityTracerName)(1, 1, 2),
            expected_boundary_tracers[0]
                .mass(gwm::state::kSpecificHumidityTracerName)(1, 1, 2),
            1.0e-5f);
  TEST_NEAR(boundary_warm_rain_tracers[0]
                .mass(gwm::state::kCloudWaterTracerName)(0, 0, 0),
            0.020f, 1.0e-6f);
  TEST_NEAR(boundary_warm_rain_tracers[0]
                .mass(gwm::state::kRainWaterTracerName)(0, 0, 0),
            0.010f, 1.0e-6f);

  qv.fill(0.0f);
  qc.fill(0.030f);
  qr.fill(0.040f);
  tracer_boundary_updater.apply(boundary_warm_rain_tracers, layout, 1800.0f);
  const auto expected_late_boundary_tracers = make_warm_rain_tracers_from_analysis(
      late_boundary_analysis, expected_late_boundary_states, layout,
      "expected_late_boundary_warm_rain");
  TEST_NEAR(boundary_warm_rain_tracers[0]
                .mass(gwm::state::kSpecificHumidityTracerName)(0, 0, 0),
            expected_late_boundary_tracers[0]
                .mass(gwm::state::kSpecificHumidityTracerName)(0, 0, 0),
            1.0e-5f);
  TEST_NEAR(boundary_warm_rain_tracers[0]
                .mass(gwm::state::kSpecificHumidityTracerName)(1, 1, 2),
            expected_late_boundary_tracers[0]
                .mass(gwm::state::kSpecificHumidityTracerName)(1, 1, 2),
            1.0e-5f);
  TEST_NEAR(boundary_warm_rain_tracers[0]
                .mass(gwm::state::kCloudWaterTracerName)(0, 0, 0),
            0.030f, 1.0e-6f);
  TEST_NEAR(boundary_warm_rain_tracers[0]
                .mass(gwm::state::kRainWaterTracerName)(0, 0, 0),
            0.040f, 1.0e-6f);

  const auto boundary_diag = diagnose_prepared_case_balance(
      late_boundary_analysis, metrics, boundary_states, layout,
      &boundary_warm_rain_tracers);
  check_balance_diag(boundary_diag);

  return 0;
}
