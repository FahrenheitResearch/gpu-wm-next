#include "gwm/ingest/prepared_case_init.hpp"

#include <algorithm>
#include <cmath>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/core/dry_thermo.hpp"
#include "gwm/dycore/boundary_conditions.hpp"
#include "gwm/dycore/passive_tracer.hpp"

namespace gwm::ingest {

namespace {

std::size_t linear_index_3d(int i, int j, int k, int nx, int ny) {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

real bundle_value_3d(const FieldBundle3D& bundle, const char* name, int i, int j,
                     int k, int nx, int ny) {
  const auto found = bundle.values.find(name);
  gwm::require(found != bundle.values.end(),
               std::string("Missing 3-D canonical field: ") + name);
  return found->second[linear_index_3d(i, j, k, nx, ny)];
}

}  // namespace

std::vector<domain::SubdomainDescriptor> make_prepared_case_layout(
    const AnalysisStateIR& analysis, const PreparedCaseInitConfig& config) {
  validate_runtime_analysis(analysis);
  gwm::require(config.halo > 0, "Prepared-case halo must be positive");
  gwm::require(config.ranks_x > 0 && config.ranks_y > 0,
               "Prepared-case rank counts must be positive");
  return comm::VirtualRankLayout::build(
      analysis.grid.nx, analysis.grid.ny, analysis.grid.nz, config.halo,
      config.ranks_x, config.ranks_y, config.periodic_x, config.periodic_y);
}

domain::GridMetrics make_prepared_case_metrics(
    const AnalysisStateIR& analysis, const PreparedCaseInitConfig& config) {
  validate_runtime_analysis(analysis);
  const auto found = analysis.static_surface.values.find("terrain_height");
  gwm::require(found != analysis.static_surface.values.end(),
               "Prepared analysis state requires terrain_height");
  return domain::GridMetrics::make_hybrid_height(
      analysis.grid.nx, analysis.grid.ny, analysis.grid.nz, analysis.grid.dx,
      analysis.grid.dy, analysis.grid.z_top, found->second,
      config.terrain_taper_eta, config.periodic_x, config.periodic_y);
}

std::vector<dycore::DryState> make_dry_states_from_analysis(
    const AnalysisStateIR& analysis,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const std::string& label_prefix) {
  validate_runtime_analysis(analysis);
  gwm::require(!layout.empty(), "Prepared-case layout must not be empty");

  const auto& pressure = analysis.atmosphere.values.at("air_pressure");
  const auto& temperature = analysis.atmosphere.values.at("air_temperature");
  const auto& u_wind = analysis.atmosphere.values.at("u_wind");
  const auto& v_wind = analysis.atmosphere.values.at("v_wind");
  const auto& w_wind = analysis.atmosphere.values.at("w_wind");

  std::vector<real> rho_global(analysis.grid.cell_count_3d(), 0.0f);
  std::vector<real> theta_global(analysis.grid.cell_count_3d(), 0.0f);
  for (int k = 0; k < analysis.grid.nz; ++k) {
    for (int j = 0; j < analysis.grid.ny; ++j) {
      for (int i = 0; i < analysis.grid.nx; ++i) {
        const auto idx =
            linear_index_3d(i, j, k, analysis.grid.nx, analysis.grid.ny);
        const real p = pressure[idx];
        const real t = temperature[idx];
        rho_global[idx] =
            p / std::max<real>(gwm::core::kDryGasConstant * t, 1.0e-3f);
        theta_global[idx] =
            t * std::pow(gwm::core::kReferencePressure /
                             std::max<real>(p, 1.0f),
                         gwm::core::kKappa);
      }
    }
  }

  std::vector<dycore::DryState> states;
  states.reserve(layout.size());
  for (const auto& desc : layout) {
    dycore::DryState state(desc.nx_local(), desc.ny_local(), desc.nz, desc.halo,
                           label_prefix + "_rank_" + std::to_string(desc.rank));
    state.fill_constant(0.0f, 300.0f, 0.0f, 0.0f, 0.0f);

    for (int j = 0; j < desc.ny_local(); ++j) {
      const int j_global = desc.j_begin + j;
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k = 0; k < desc.nz; ++k) {
          const auto idx =
              linear_index_3d(i_global, j_global, k, analysis.grid.nx,
                              analysis.grid.ny);
          state.rho_d(i, j, k) = rho_global[idx];
          state.rho_theta_m(i, j, k) = rho_global[idx] * theta_global[idx];
        }
      }
    }

    for (int j = 0; j < desc.ny_local(); ++j) {
      const int j_global = desc.j_begin + j;
      for (int i_face = 0; i_face < desc.nx_local() + 1; ++i_face) {
        const int i_right =
            std::clamp(desc.i_begin + i_face, 0, analysis.grid.nx - 1);
        const int i_left = std::clamp(i_right - 1, 0, analysis.grid.nx - 1);
        for (int k = 0; k < desc.nz; ++k) {
          const auto idx_left =
              linear_index_3d(i_left, j_global, k, analysis.grid.nx,
                              analysis.grid.ny);
          const auto idx_right =
              linear_index_3d(i_right, j_global, k, analysis.grid.nx,
                              analysis.grid.ny);
          const real rho_face = 0.5f * (rho_global[idx_left] + rho_global[idx_right]);
          const real u_face = 0.5f * (u_wind[idx_left] + u_wind[idx_right]);
          state.mom_u.storage()(i_face, j, k) = rho_face * u_face;
        }
      }
    }

    for (int j_face = 0; j_face < desc.ny_local() + 1; ++j_face) {
      const int j_north =
          std::clamp(desc.j_begin + j_face, 0, analysis.grid.ny - 1);
      const int j_south = std::clamp(j_north - 1, 0, analysis.grid.ny - 1);
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k = 0; k < desc.nz; ++k) {
          const auto idx_south =
              linear_index_3d(i_global, j_south, k, analysis.grid.nx,
                              analysis.grid.ny);
          const auto idx_north =
              linear_index_3d(i_global, j_north, k, analysis.grid.nx,
                              analysis.grid.ny);
          const real rho_face = 0.5f * (rho_global[idx_south] + rho_global[idx_north]);
          const real v_face = 0.5f * (v_wind[idx_south] + v_wind[idx_north]);
          state.mom_v.storage()(i, j_face, k) = rho_face * v_face;
        }
      }
    }

    for (int j = 0; j < desc.ny_local(); ++j) {
      const int j_global = desc.j_begin + j;
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k_face = 0; k_face < desc.nz + 1; ++k_face) {
          const int k_upper = std::clamp(k_face, 0, analysis.grid.nz - 1);
          const int k_lower = std::clamp(k_upper - 1, 0, analysis.grid.nz - 1);
          const auto idx_lower =
              linear_index_3d(i_global, j_global, k_lower, analysis.grid.nx,
                              analysis.grid.ny);
          const auto idx_upper =
              linear_index_3d(i_global, j_global, k_upper, analysis.grid.nx,
                              analysis.grid.ny);
          const real rho_face = 0.5f * (rho_global[idx_lower] + rho_global[idx_upper]);
          const real w_face = 0.5f * (w_wind[idx_lower] + w_wind[idx_upper]);
          state.mom_w.storage()(i, j, k_face) = rho_face * w_face;
        }
      }
    }

    states.push_back(std::move(state));
  }

  return states;
}

std::vector<state::TracerState> make_specific_humidity_tracers_from_analysis(
    const AnalysisStateIR& analysis,
    const std::vector<dycore::DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const std::string& label_prefix) {
  validate_runtime_analysis(analysis);
  const auto& specific_humidity =
      analysis.atmosphere.values.at("specific_humidity");
  return dycore::make_specific_humidity_tracers_from_global_field(
      dry_states, layout, analysis.grid.nx, analysis.grid.ny, specific_humidity,
      label_prefix);
}

PreparedCaseBoundaryUpdater::PreparedCaseBoundaryUpdater(
    AnalysisStateIR analysis, BoundaryCacheIR cache, PreparedCaseInitConfig config)
    : analysis_(std::move(analysis)),
      cache_(std::move(cache)),
      config_(config) {}

void PreparedCaseBoundaryUpdater::set_step_start_time(real step_start_time_seconds) {
  step_start_time_seconds_ = std::max(step_start_time_seconds, 0.0f);
}

void PreparedCaseBoundaryUpdater::apply(
    std::vector<dycore::DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout, real sim_time) {
  if (cache_.empty()) {
    return;
  }

  const auto snapshot = interpolate_boundary_snapshot(
      cache_, static_cast<int>(std::lround(step_start_time_seconds_ + sim_time)));
  AnalysisStateIR boundary_analysis = analysis_;
  boundary_analysis.valid_time_utc = snapshot.valid_time_utc;
  boundary_analysis.forecast_offset_seconds = snapshot.forecast_offset_seconds;
  boundary_analysis.atmosphere = snapshot.atmosphere;
  boundary_analysis.surface = snapshot.surface;
  const auto boundary_states =
      make_dry_states_from_analysis(boundary_analysis, layout, "boundary_ref");
  dycore::apply_reference_boundaries(states, boundary_states, layout);
}

PreparedCaseTracerBoundaryUpdater::PreparedCaseTracerBoundaryUpdater(
    AnalysisStateIR analysis, BoundaryCacheIR cache, PreparedCaseInitConfig config)
    : analysis_(std::move(analysis)),
      cache_(std::move(cache)),
      config_(config) {}

void PreparedCaseTracerBoundaryUpdater::set_step_start_time(
    real step_start_time_seconds) {
  step_start_time_seconds_ = std::max(step_start_time_seconds, 0.0f);
}

void PreparedCaseTracerBoundaryUpdater::apply(
    std::vector<state::TracerState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout, real sim_time) {
  if (cache_.empty()) {
    return;
  }

  const auto snapshot = interpolate_boundary_snapshot(
      cache_, static_cast<int>(std::lround(step_start_time_seconds_ + sim_time)));
  AnalysisStateIR boundary_analysis = analysis_;
  boundary_analysis.valid_time_utc = snapshot.valid_time_utc;
  boundary_analysis.forecast_offset_seconds = snapshot.forecast_offset_seconds;
  boundary_analysis.atmosphere = snapshot.atmosphere;
  boundary_analysis.surface = snapshot.surface;
  const auto boundary_dry_states =
      make_dry_states_from_analysis(boundary_analysis, layout, "boundary_ref");
  const auto boundary_tracers = make_specific_humidity_tracers_from_analysis(
      boundary_analysis, boundary_dry_states, layout, "boundary_qv");
  dycore::apply_reference_boundaries(states, boundary_tracers, layout);
}

}  // namespace gwm::ingest
