#include "gwm/ingest/prepared_case_init.hpp"

#include <algorithm>
#include <cmath>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/core/dry_thermo.hpp"
#include "gwm/dycore/boundary_conditions.hpp"
#include "gwm/dycore/passive_tracer.hpp"

namespace gwm::ingest {

namespace {

constexpr real kPreparedCaseGravity = 9.81f;

std::size_t linear_index_3d(int i, int j, int k, int nx, int ny) {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

void require_metrics_compatible(const AnalysisStateIR& analysis,
                                const domain::GridMetrics& metrics) {
  gwm::require(metrics.nx == analysis.grid.nx && metrics.ny == analysis.grid.ny &&
                   metrics.nz == analysis.grid.nz,
               "Prepared-case metrics must match analysis grid dimensions");
}

int wrap_or_clamp_index(int idx, int count, bool periodic) {
  if (!periodic) {
    return std::clamp(idx, 0, count - 1);
  }
  const int mod = idx % count;
  return mod < 0 ? mod + count : mod;
}

real rho_theta_m_from_pressure(real pressure) {
  real rho_theta =
      (gwm::core::kReferencePressure / gwm::core::kDryGasConstant) *
      std::pow(
          std::max<real>(pressure / gwm::core::kReferencePressure, 1.0e-6f),
          1.0f - gwm::core::kKappa);
  for (int iter = 0; iter < 2; ++iter) {
    const real pressure_eval =
        gwm::core::dry_pressure_from_rho_theta_m(1.0f, rho_theta);
    const real ratio =
        std::max<real>(pressure / std::max<real>(pressure_eval, 1.0f), 1.0e-6f);
    rho_theta *= std::pow(ratio, 1.0f - gwm::core::kKappa);
  }
  return rho_theta;
}

struct BalancedPreparedCaseColumns {
  std::vector<real> pressure;
  std::vector<real> theta;
  std::vector<real> rho_d;
  std::vector<real> rho_theta_m;
};

std::size_t linear_index_xface(int i_face, int j, int k, int nx, int ny) {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(nx + 1) +
         static_cast<std::size_t>(i_face);
}

std::size_t linear_index_yface(int i, int j_face, int k, int nx, int ny) {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny + 1) +
          static_cast<std::size_t>(j_face)) *
             static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

std::size_t linear_index_zface(int i, int j, int k_face, int nx, int ny) {
  return (static_cast<std::size_t>(k_face) * static_cast<std::size_t>(ny) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

std::size_t linear_index_2d(int i, int j, int nx) {
  return static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

struct PreparedCaseFaceMomenta {
  std::vector<real> mom_u;
  std::vector<real> mom_v;
  std::vector<real> mom_w;
};

std::vector<real> compute_column_mean_divergence(
    const std::vector<real>& mom_u, const std::vector<real>& mom_v,
    const AnalysisStateIR& analysis, const domain::GridMetrics& metrics) {
  const int nx = analysis.grid.nx;
  const int ny = analysis.grid.ny;
  const int nz = analysis.grid.nz;
  std::vector<real> column_mean_divergence(
      static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny), 0.0f);
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      real column_mass_bias = 0.0f;
      real column_height = 0.0f;
      for (int k = 0; k < nz; ++k) {
        const real hdiv =
            (mom_u[linear_index_xface(i + 1, j, k, nx, ny)] -
             mom_u[linear_index_xface(i, j, k, nx, ny)]) /
                static_cast<real>(analysis.grid.dx) +
            (mom_v[linear_index_yface(i, j + 1, k, nx, ny)] -
             mom_v[linear_index_yface(i, j, k, nx, ny)]) /
                static_cast<real>(analysis.grid.dy);
        const real inv_dz = metrics.inv_dz_cell(i, j, k);
        const real dz = 1.0f / std::max(inv_dz, 1.0e-6f);
        column_mass_bias += hdiv * dz;
        column_height += dz;
      }
      column_mean_divergence[linear_index_2d(i, j, nx)] =
          column_mass_bias / std::max(column_height, 1.0e-6f);
    }
  }
  return column_mean_divergence;
}

std::vector<real> solve_barotropic_projection_potential(
    const std::vector<real>& rhs, int nx, int ny, real dx, real dy) {
  std::vector<real> phi(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny),
                        0.0f);
  if (nx == 0 || ny == 0) {
    return phi;
  }

  const real inv_dx2 = 1.0f / std::max(dx * dx, 1.0e-6f);
  const real inv_dy2 = 1.0f / std::max(dy * dy, 1.0e-6f);
  const real denom = 2.0f * (inv_dx2 + inv_dy2);
  if (denom <= 0.0f) {
    return phi;
  }

  constexpr int kMaxIterations = 2000;
  constexpr real kOmega = 1.85f;
  constexpr real kTolerance = 1.0e-10f;
  for (int iter = 0; iter < kMaxIterations; ++iter) {
    real max_update = 0.0f;
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        const real phi_west = i > 0 ? phi[linear_index_2d(i - 1, j, nx)] : 0.0f;
        const real phi_east =
            i + 1 < nx ? phi[linear_index_2d(i + 1, j, nx)] : 0.0f;
        const real phi_south = j > 0 ? phi[linear_index_2d(i, j - 1, nx)] : 0.0f;
        const real phi_north =
            j + 1 < ny ? phi[linear_index_2d(i, j + 1, nx)] : 0.0f;
        const auto idx = linear_index_2d(i, j, nx);
        const real phi_star =
            ((phi_west + phi_east) * inv_dx2 +
             (phi_south + phi_north) * inv_dy2 -
             rhs[idx]) /
            denom;
        const real updated = phi[idx] + kOmega * (phi_star - phi[idx]);
        max_update = std::max(max_update, std::fabs(updated - phi[idx]));
        phi[idx] = updated;
      }
    }
    if (max_update <= kTolerance) {
      break;
    }
  }
  return phi;
}

PreparedCaseFaceMomenta build_projected_prepared_case_face_momenta(
    const AnalysisStateIR& analysis, const domain::GridMetrics& metrics,
    const BalancedPreparedCaseColumns& balanced) {
  const int nx = analysis.grid.nx;
  const int ny = analysis.grid.ny;
  const int nz = analysis.grid.nz;
  const auto& u_wind = analysis.atmosphere.values.at("u_wind");
  const auto& v_wind = analysis.atmosphere.values.at("v_wind");

  PreparedCaseFaceMomenta face_momenta{};
  face_momenta.mom_u.resize(static_cast<std::size_t>(nx + 1) *
                                static_cast<std::size_t>(ny) *
                                static_cast<std::size_t>(nz),
                            0.0f);
  face_momenta.mom_v.resize(static_cast<std::size_t>(nx) *
                                static_cast<std::size_t>(ny + 1) *
                                static_cast<std::size_t>(nz),
                            0.0f);
  face_momenta.mom_w.resize(static_cast<std::size_t>(nx) *
                                static_cast<std::size_t>(ny) *
                                static_cast<std::size_t>(nz + 1),
                            0.0f);

  for (int j = 0; j < ny; ++j) {
    for (int i_face = 0; i_face < nx + 1; ++i_face) {
      const int i_right = wrap_or_clamp_index(i_face, nx, metrics.periodic_x);
      const int i_left = wrap_or_clamp_index(i_face - 1, nx, metrics.periodic_x);
      for (int k = 0; k < nz; ++k) {
        const auto idx_left = linear_index_3d(i_left, j, k, nx, ny);
        const auto idx_right = linear_index_3d(i_right, j, k, nx, ny);
        const real rho_face =
            0.5f * (balanced.rho_d[idx_left] + balanced.rho_d[idx_right]);
        const real u_face = 0.5f * (u_wind[idx_left] + u_wind[idx_right]);
        face_momenta.mom_u[linear_index_xface(i_face, j, k, nx, ny)] =
            rho_face * u_face;
      }
    }
  }

  for (int j_face = 0; j_face < ny + 1; ++j_face) {
    const int j_north = wrap_or_clamp_index(j_face, ny, metrics.periodic_y);
    const int j_south = wrap_or_clamp_index(j_face - 1, ny, metrics.periodic_y);
    for (int i = 0; i < nx; ++i) {
      for (int k = 0; k < nz; ++k) {
        const auto idx_south = linear_index_3d(i, j_south, k, nx, ny);
        const auto idx_north = linear_index_3d(i, j_north, k, nx, ny);
        const real rho_face =
            0.5f * (balanced.rho_d[idx_south] + balanced.rho_d[idx_north]);
        const real v_face = 0.5f * (v_wind[idx_south] + v_wind[idx_north]);
        face_momenta.mom_v[linear_index_yface(i, j_face, k, nx, ny)] =
            rho_face * v_face;
      }
    }
  }

  if (!metrics.periodic_x && !metrics.periodic_y) {
    constexpr int kProjectionRefinements = 3;
    constexpr real kProjectionTolerance = 1.0e-8f;
    for (int refinement = 0; refinement < kProjectionRefinements; ++refinement) {
      const auto column_mean_divergence = compute_column_mean_divergence(
          face_momenta.mom_u, face_momenta.mom_v, analysis, metrics);
      real max_abs_mean_divergence = 0.0f;
      for (const real value : column_mean_divergence) {
        max_abs_mean_divergence =
            std::max(max_abs_mean_divergence, std::fabs(value));
      }
      if (max_abs_mean_divergence <= kProjectionTolerance) {
        break;
      }

      const auto phi = solve_barotropic_projection_potential(
          column_mean_divergence, nx, ny, static_cast<real>(analysis.grid.dx),
          static_cast<real>(analysis.grid.dy));
      for (int j = 0; j < ny; ++j) {
        for (int i_face = 0; i_face < nx + 1; ++i_face) {
          const real phi_left =
              i_face > 0 ? phi[linear_index_2d(i_face - 1, j, nx)] : 0.0f;
          const real phi_right =
              i_face < nx ? phi[linear_index_2d(i_face, j, nx)] : 0.0f;
          const real correction =
              (phi_right - phi_left) / static_cast<real>(analysis.grid.dx);
          for (int k = 0; k < nz; ++k) {
            face_momenta.mom_u[linear_index_xface(i_face, j, k, nx, ny)] -=
                correction;
          }
        }
      }

      for (int j_face = 0; j_face < ny + 1; ++j_face) {
        for (int i = 0; i < nx; ++i) {
          const real phi_south =
              j_face > 0 ? phi[linear_index_2d(i, j_face - 1, nx)] : 0.0f;
          const real phi_north =
              j_face < ny ? phi[linear_index_2d(i, j_face, nx)] : 0.0f;
          const real correction =
              (phi_north - phi_south) / static_cast<real>(analysis.grid.dy);
          for (int k = 0; k < nz; ++k) {
            face_momenta.mom_v[linear_index_yface(i, j_face, k, nx, ny)] -=
                correction;
          }
        }
      }
    }
  }

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      face_momenta.mom_w[linear_index_zface(i, j, 0, nx, ny)] = 0.0f;

      real column_mass_bias = 0.0f;
      real column_height = 0.0f;
      for (int k = 0; k < nz; ++k) {
        const real hdiv =
            (face_momenta.mom_u[linear_index_xface(i + 1, j, k, nx, ny)] -
             face_momenta.mom_u[linear_index_xface(i, j, k, nx, ny)]) /
                static_cast<real>(analysis.grid.dx) +
            (face_momenta.mom_v[linear_index_yface(i, j + 1, k, nx, ny)] -
             face_momenta.mom_v[linear_index_yface(i, j, k, nx, ny)]) /
                static_cast<real>(analysis.grid.dy);
        const real inv_dz = metrics.inv_dz_cell(i, j, k);
        const real dz = 1.0f / std::max(inv_dz, 1.0e-6f);
        column_mass_bias += hdiv * dz;
        column_height += dz;
      }

      const real mean_divergence =
          column_mass_bias / std::max(column_height, 1.0e-6f);
      for (int k = 0; k < nz; ++k) {
        const real hdiv =
            (face_momenta.mom_u[linear_index_xface(i + 1, j, k, nx, ny)] -
             face_momenta.mom_u[linear_index_xface(i, j, k, nx, ny)]) /
                static_cast<real>(analysis.grid.dx) +
            (face_momenta.mom_v[linear_index_yface(i, j + 1, k, nx, ny)] -
             face_momenta.mom_v[linear_index_yface(i, j, k, nx, ny)]) /
                static_cast<real>(analysis.grid.dy);
        const real inv_dz = metrics.inv_dz_cell(i, j, k);
        const real dz = 1.0f / std::max(inv_dz, 1.0e-6f);
        face_momenta.mom_w[linear_index_zface(i, j, k + 1, nx, ny)] =
            face_momenta.mom_w[linear_index_zface(i, j, k, nx, ny)] -
            (hdiv - mean_divergence) * dz;
      }
      face_momenta.mom_w[linear_index_zface(i, j, nz, nx, ny)] = 0.0f;
    }
  }

  return face_momenta;
}

BalancedPreparedCaseColumns build_balanced_prepared_case_columns(
    const AnalysisStateIR& analysis, const domain::GridMetrics& metrics) {
  require_metrics_compatible(analysis, metrics);

  const int nx = analysis.grid.nx;
  const int ny = analysis.grid.ny;
  const int nz = analysis.grid.nz;
  const auto cell_count = analysis.grid.cell_count_3d();
  const auto& pressure_src = analysis.atmosphere.values.at("air_pressure");
  const auto& temperature_src = analysis.atmosphere.values.at("air_temperature");

  BalancedPreparedCaseColumns columns{};
  columns.pressure.resize(cell_count, 0.0f);
  columns.theta.resize(cell_count, 0.0f);
  columns.rho_d.resize(cell_count, 0.0f);
  columns.rho_theta_m.resize(cell_count, 0.0f);

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const auto idx0 = linear_index_3d(i, j, 0, nx, ny);
      const real p0 = std::max<real>(pressure_src[idx0], 1.0f);
      const real t0 = temperature_src[idx0];
      const real theta0 =
          t0 * std::pow(gwm::core::kReferencePressure / p0,
                        gwm::core::kKappa);
      columns.theta[idx0] = theta0;
      columns.pressure[idx0] = p0;
      columns.rho_d[idx0] = gwm::core::dry_rho_from_pressure_theta(p0, theta0);
      columns.rho_theta_m[idx0] = rho_theta_m_from_pressure(p0);

      for (int k = 1; k < nz; ++k) {
        const auto idx = linear_index_3d(i, j, k, nx, ny);
        const auto idx_prev = linear_index_3d(i, j, k - 1, nx, ny);
        const real p_src = std::max<real>(pressure_src[idx], 1.0f);
        const real temperature_k = temperature_src[idx];
        const real theta_k =
            temperature_k *
            std::pow(gwm::core::kReferencePressure / p_src,
                     gwm::core::kKappa);
        columns.theta[idx] = theta_k;

        const real dz =
            metrics.z_center(i, j, k) - metrics.z_center(i, j, k - 1);
        real p_next = std::max<real>(
            1.0f, columns.pressure[idx_prev] -
                      kPreparedCaseGravity * columns.rho_d[idx_prev] * dz);
        for (int iter = 0; iter < 4; ++iter) {
          const real rho_next =
              gwm::core::dry_rho_from_pressure_theta(p_next, theta_k);
          p_next = columns.pressure[idx_prev] -
                   kPreparedCaseGravity *
                       0.5f * (columns.rho_d[idx_prev] + rho_next) * dz;
          p_next = std::max<real>(p_next, 1.0f);
        }

        columns.pressure[idx] = p_next;
        columns.rho_d[idx] =
            gwm::core::dry_rho_from_pressure_theta(p_next, theta_k);
        columns.rho_theta_m[idx] = rho_theta_m_from_pressure(p_next);
      }
    }
  }

  return columns;
}

AnalysisStateIR make_boundary_analysis_from_snapshot(
    const AnalysisStateIR& analysis, const BoundarySnapshotIR& snapshot) {
  AnalysisStateIR boundary_analysis = analysis;
  boundary_analysis.valid_time_utc = snapshot.valid_time_utc;
  boundary_analysis.forecast_offset_seconds = snapshot.forecast_offset_seconds;
  boundary_analysis.atmosphere = snapshot.atmosphere;
  boundary_analysis.surface = snapshot.surface;
  return boundary_analysis;
}

void refresh_boundary_projection_cache(
    const AnalysisStateIR& analysis, const BoundaryCacheIR& cache,
    const domain::GridMetrics& metrics, real step_start_time_seconds,
    real sim_time, int& cached_forecast_offset_seconds,
    BoundarySnapshotIR& cached_snapshot, std::vector<real>& cached_rho_d,
    std::vector<real>& cached_rho_theta_m, std::vector<real>& cached_mom_u,
    std::vector<real>& cached_mom_v, std::vector<real>& cached_mom_w) {
  const int forecast_offset_seconds =
      static_cast<int>(std::lround(step_start_time_seconds + sim_time));
  if (forecast_offset_seconds == cached_forecast_offset_seconds) {
    return;
  }

  cached_snapshot = interpolate_boundary_snapshot(cache, forecast_offset_seconds);
  const auto boundary_analysis =
      make_boundary_analysis_from_snapshot(analysis, cached_snapshot);
  auto balanced = build_balanced_prepared_case_columns(boundary_analysis, metrics);
  auto face_momenta = build_projected_prepared_case_face_momenta(
      boundary_analysis, metrics, balanced);

  cached_forecast_offset_seconds = forecast_offset_seconds;
  cached_rho_d = std::move(balanced.rho_d);
  cached_rho_theta_m = std::move(balanced.rho_theta_m);
  cached_mom_u = std::move(face_momenta.mom_u);
  cached_mom_v = std::move(face_momenta.mom_v);
  cached_mom_w = std::move(face_momenta.mom_w);
}

void apply_cached_dry_boundary_strips(
    std::vector<dycore::DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout, int global_nx,
    int global_ny, const std::vector<real>& rho_d,
    const std::vector<real>& rho_theta_m, const std::vector<real>& mom_u,
    const std::vector<real>& mom_v, const std::vector<real>& mom_w) {
  gwm::require(states.size() == layout.size(),
               "State/layout size mismatch in apply_cached_dry_boundary_strips");

  for (std::size_t rank = 0; rank < layout.size(); ++rank) {
    auto& state = states[rank];
    const auto& desc = layout[rank];

    if (dycore::touches_west_boundary(desc)) {
      const int i_global = desc.i_begin;
      for (int j = 0; j < desc.ny_local(); ++j) {
        const int j_global = desc.j_begin + j;
        for (int k = 0; k < desc.nz; ++k) {
          const auto idx = linear_index_3d(i_global, j_global, k, global_nx,
                                           global_ny);
          state.rho_d(0, j, k) = rho_d[idx];
          state.rho_theta_m(0, j, k) = rho_theta_m[idx];
          state.mom_u.storage()(0, j, k) =
              mom_u[linear_index_xface(i_global, j_global, k, global_nx,
                                       global_ny)];
        }
      }
      for (int j_face = 0; j_face < desc.ny_local() + 1; ++j_face) {
        const int j_global = desc.j_begin + j_face;
        for (int k = 0; k < desc.nz; ++k) {
          state.mom_v.storage()(0, j_face, k) =
              mom_v[linear_index_yface(i_global, j_global, k, global_nx,
                                       global_ny)];
        }
      }
      for (int j = 0; j < desc.ny_local(); ++j) {
        const int j_global = desc.j_begin + j;
        for (int k_face = 0; k_face <= desc.nz; ++k_face) {
          state.mom_w.storage()(0, j, k_face) =
              mom_w[linear_index_zface(i_global, j_global, k_face, global_nx,
                                       global_ny)];
        }
      }
    }

    if (dycore::touches_east_boundary(desc)) {
      const int i_local = desc.nx_local() - 1;
      const int i_global = desc.i_begin + i_local;
      for (int j = 0; j < desc.ny_local(); ++j) {
        const int j_global = desc.j_begin + j;
        for (int k = 0; k < desc.nz; ++k) {
          const auto idx = linear_index_3d(i_global, j_global, k, global_nx,
                                           global_ny);
          state.rho_d(i_local, j, k) = rho_d[idx];
          state.rho_theta_m(i_local, j, k) = rho_theta_m[idx];
          state.mom_u.storage()(desc.nx_local(), j, k) =
              mom_u[linear_index_xface(desc.i_begin + desc.nx_local(), j_global,
                                       k, global_nx, global_ny)];
        }
      }
      for (int j_face = 0; j_face < desc.ny_local() + 1; ++j_face) {
        const int j_global = desc.j_begin + j_face;
        for (int k = 0; k < desc.nz; ++k) {
          state.mom_v.storage()(i_local, j_face, k) =
              mom_v[linear_index_yface(i_global, j_global, k, global_nx,
                                       global_ny)];
        }
      }
      for (int j = 0; j < desc.ny_local(); ++j) {
        const int j_global = desc.j_begin + j;
        for (int k_face = 0; k_face <= desc.nz; ++k_face) {
          state.mom_w.storage()(i_local, j, k_face) =
              mom_w[linear_index_zface(i_global, j_global, k_face, global_nx,
                                       global_ny)];
        }
      }
    }

    if (dycore::touches_south_boundary(desc)) {
      const int j_global = desc.j_begin;
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k = 0; k < desc.nz; ++k) {
          const auto idx = linear_index_3d(i_global, j_global, k, global_nx,
                                           global_ny);
          state.rho_d(i, 0, k) = rho_d[idx];
          state.rho_theta_m(i, 0, k) = rho_theta_m[idx];
          state.mom_v.storage()(i, 0, k) =
              mom_v[linear_index_yface(i_global, j_global, k, global_nx,
                                       global_ny)];
        }
      }
      for (int i_face = 0; i_face < desc.nx_local() + 1; ++i_face) {
        const int i_global = desc.i_begin + i_face;
        for (int k = 0; k < desc.nz; ++k) {
          state.mom_u.storage()(i_face, 0, k) =
              mom_u[linear_index_xface(i_global, j_global, k, global_nx,
                                       global_ny)];
        }
      }
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k_face = 0; k_face <= desc.nz; ++k_face) {
          state.mom_w.storage()(i, 0, k_face) =
              mom_w[linear_index_zface(i_global, j_global, k_face, global_nx,
                                       global_ny)];
        }
      }
    }

    if (dycore::touches_north_boundary(desc)) {
      const int j_local = desc.ny_local() - 1;
      const int j_global = desc.j_begin + j_local;
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k = 0; k < desc.nz; ++k) {
          const auto idx = linear_index_3d(i_global, j_global, k, global_nx,
                                           global_ny);
          state.rho_d(i, j_local, k) = rho_d[idx];
          state.rho_theta_m(i, j_local, k) = rho_theta_m[idx];
          state.mom_v.storage()(i, desc.ny_local(), k) =
              mom_v[linear_index_yface(i_global, desc.j_begin + desc.ny_local(),
                                       k, global_nx, global_ny)];
        }
      }
      for (int i_face = 0; i_face < desc.nx_local() + 1; ++i_face) {
        const int i_global = desc.i_begin + i_face;
        for (int k = 0; k < desc.nz; ++k) {
          state.mom_u.storage()(i_face, j_local, k) =
              mom_u[linear_index_xface(i_global, j_global, k, global_nx,
                                       global_ny)];
        }
      }
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k_face = 0; k_face <= desc.nz; ++k_face) {
          state.mom_w.storage()(i, j_local, k_face) =
              mom_w[linear_index_zface(i_global, j_global, k_face, global_nx,
                                       global_ny)];
        }
      }
    }
  }
}

void apply_cached_tracer_boundary_strips(
    std::vector<state::TracerState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const BoundarySnapshotIR& snapshot, int global_nx, int global_ny,
    const std::vector<real>& rho_d) {
  gwm::require(
      states.size() == layout.size(),
      "Tracer/layout size mismatch in apply_cached_tracer_boundary_strips");

  for (std::size_t rank = 0; rank < layout.size(); ++rank) {
    auto& state = states[rank];
    const auto& desc = layout[rank];
    for (std::size_t tracer_index = 0; tracer_index < state.size();
         ++tracer_index) {
      const auto& spec = state.registry().at(static_cast<int>(tracer_index));
      const auto field_it = snapshot.atmosphere.values.find(spec.name);
      if (field_it == snapshot.atmosphere.values.end()) {
        continue;
      }

      const auto& mixing_ratio = field_it->second;
      auto& rho_q = state.mass(static_cast<int>(tracer_index));
      if (dycore::touches_west_boundary(desc)) {
        const int i_global = desc.i_begin;
        for (int j = 0; j < desc.ny_local(); ++j) {
          const int j_global = desc.j_begin + j;
          for (int k = 0; k < desc.nz; ++k) {
            const auto idx = linear_index_3d(i_global, j_global, k, global_nx,
                                             global_ny);
            rho_q(0, j, k) = rho_d[idx] * mixing_ratio[idx];
          }
        }
      }

      if (dycore::touches_east_boundary(desc)) {
        const int i_local = desc.nx_local() - 1;
        const int i_global = desc.i_begin + i_local;
        for (int j = 0; j < desc.ny_local(); ++j) {
          const int j_global = desc.j_begin + j;
          for (int k = 0; k < desc.nz; ++k) {
            const auto idx = linear_index_3d(i_global, j_global, k, global_nx,
                                             global_ny);
            rho_q(i_local, j, k) = rho_d[idx] * mixing_ratio[idx];
          }
        }
      }

      if (dycore::touches_south_boundary(desc)) {
        const int j_global = desc.j_begin;
        for (int i = 0; i < desc.nx_local(); ++i) {
          const int i_global = desc.i_begin + i;
          for (int k = 0; k < desc.nz; ++k) {
            const auto idx = linear_index_3d(i_global, j_global, k, global_nx,
                                             global_ny);
            rho_q(i, 0, k) = rho_d[idx] * mixing_ratio[idx];
          }
        }
      }

      if (dycore::touches_north_boundary(desc)) {
        const int j_local = desc.ny_local() - 1;
        const int j_global = desc.j_begin + j_local;
        for (int i = 0; i < desc.nx_local(); ++i) {
          const int i_global = desc.i_begin + i;
          for (int k = 0; k < desc.nz; ++k) {
            const auto idx = linear_index_3d(i_global, j_global, k, global_nx,
                                             global_ny);
            rho_q(i, j_local, k) = rho_d[idx] * mixing_ratio[idx];
          }
        }
      }
    }
  }
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
    const AnalysisStateIR& analysis, const domain::GridMetrics& metrics,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const std::string& label_prefix) {
  validate_runtime_analysis(analysis);
  gwm::require(!layout.empty(), "Prepared-case layout must not be empty");
  require_metrics_compatible(analysis, metrics);

  const auto& u_wind = analysis.atmosphere.values.at("u_wind");
  const auto& v_wind = analysis.atmosphere.values.at("v_wind");
  const auto balanced = build_balanced_prepared_case_columns(analysis, metrics);
  const auto face_momenta =
      build_projected_prepared_case_face_momenta(analysis, metrics, balanced);

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
          state.rho_d(i, j, k) = balanced.rho_d[idx];
          state.rho_theta_m(i, j, k) = balanced.rho_theta_m[idx];
        }
      }
    }

    for (int j = 0; j < desc.ny_local(); ++j) {
      const int j_global = desc.j_begin + j;
      for (int i_face = 0; i_face < desc.nx_local() + 1; ++i_face) {
        for (int k = 0; k < desc.nz; ++k) {
          state.mom_u.storage()(i_face, j, k) =
              face_momenta.mom_u[linear_index_xface(desc.i_begin + i_face, j_global,
                                                    k, analysis.grid.nx,
                                                    analysis.grid.ny)];
        }
      }
    }

    for (int j_face = 0; j_face < desc.ny_local() + 1; ++j_face) {
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k = 0; k < desc.nz; ++k) {
          state.mom_v.storage()(i, j_face, k) =
              face_momenta.mom_v[linear_index_yface(i_global, desc.j_begin + j_face,
                                                    k, analysis.grid.nx,
                                                    analysis.grid.ny)];
        }
      }
    }

    for (int j = 0; j < desc.ny_local(); ++j) {
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k_face = 0; k_face <= desc.nz; ++k_face) {
          state.mom_w.storage()(i, j, k_face) =
              face_momenta.mom_w[linear_index_zface(i_global, desc.j_begin + j,
                                                    k_face, analysis.grid.nx,
                                                    analysis.grid.ny)];
        }
      }
    }

    states.push_back(std::move(state));
  }

  return states;
}

PreparedCaseBalanceDiagnostics diagnose_prepared_case_balance(
    const AnalysisStateIR& analysis, const domain::GridMetrics& metrics,
    const std::vector<dycore::DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const std::vector<state::TracerState>* tracers) {
  validate_runtime_analysis(analysis);
  require_metrics_compatible(analysis, metrics);
  gwm::require(dry_states.size() == layout.size(),
               "Dry-state/layout size mismatch in diagnose_prepared_case_balance");
  if (tracers != nullptr) {
    gwm::require(tracers->size() == layout.size(),
                 "Tracer/layout size mismatch in diagnose_prepared_case_balance");
  }

  const auto balanced = build_balanced_prepared_case_columns(analysis, metrics);
  PreparedCaseBalanceDiagnostics diagnostics{};

  std::vector<real> rho_global(analysis.grid.cell_count_3d(), 0.0f);
  std::vector<real> theta_global(analysis.grid.cell_count_3d(), 0.0f);
  std::vector<real> pressure_global(analysis.grid.cell_count_3d(), 0.0f);

  for (std::size_t rank = 0; rank < layout.size(); ++rank) {
    const auto& desc = layout[rank];
    const auto& state = dry_states[rank];
    for (int j = 0; j < desc.ny_local(); ++j) {
      const int j_global = desc.j_begin + j;
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k = 0; k < desc.nz; ++k) {
          const auto idx =
              linear_index_3d(i_global, j_global, k, analysis.grid.nx,
                              analysis.grid.ny);
          rho_global[idx] = state.rho_d(i, j, k);
          theta_global[idx] = gwm::core::dry_theta_from_conserved(
              state.rho_d(i, j, k), state.rho_theta_m(i, j, k));
          pressure_global[idx] = gwm::core::dry_pressure_from_rho_theta_m(
              state.rho_d(i, j, k), state.rho_theta_m(i, j, k));

          const real div_x =
              (state.mom_u.storage()(i + 1, j, k) - state.mom_u.storage()(i, j, k)) /
              static_cast<real>(analysis.grid.dx);
          const real div_y =
              (state.mom_v.storage()(i, j + 1, k) - state.mom_v.storage()(i, j, k)) /
              static_cast<real>(analysis.grid.dy);
          const real div_z =
              (state.mom_w.storage()(i, j, k + 1) - state.mom_w.storage()(i, j, k)) *
              metrics.inv_dz_cell(i_global, j_global, k);
          const real div_mass = div_x + div_y + div_z;
          const real div_scale = std::max(
              {std::fabs(div_x), std::fabs(div_y), std::fabs(div_z), 1.0e-12f});
          diagnostics.max_abs_mass_divergence = std::max(
              diagnostics.max_abs_mass_divergence, std::fabs(div_mass));
          diagnostics.max_rel_mass_divergence = std::max(
              diagnostics.max_rel_mass_divergence,
              std::fabs(div_mass) / div_scale);
        }

        diagnostics.max_abs_mom_w_bottom =
            std::max(diagnostics.max_abs_mom_w_bottom,
                     std::fabs(state.mom_w.storage()(i, j, 0)));
        diagnostics.max_abs_mom_w_top =
            std::max(diagnostics.max_abs_mom_w_top,
                     std::fabs(state.mom_w.storage()(i, j, desc.nz)));
      }
    }
  }

  const auto geopotential_height =
      analysis.atmosphere.values.find("geopotential_height");
  for (int j = 0; j < analysis.grid.ny; ++j) {
    for (int i = 0; i < analysis.grid.nx; ++i) {
      for (int k = 0; k < analysis.grid.nz; ++k) {
        const auto idx =
            linear_index_3d(i, j, k, analysis.grid.nx, analysis.grid.ny);
        const real rho_expected = gwm::core::dry_rho_from_pressure_theta(
            balanced.pressure[idx], theta_global[idx]);
        diagnostics.max_rel_eos = std::max(
            diagnostics.max_rel_eos,
            std::fabs(rho_global[idx] - rho_expected) /
                std::max(std::fabs(rho_expected), 1.0f));
        diagnostics.max_rel_hydrostatic = std::max(
            diagnostics.max_rel_hydrostatic,
            std::fabs(pressure_global[idx] - balanced.pressure[idx]) /
                std::max(std::fabs(balanced.pressure[idx]), 1.0f));

        if (geopotential_height != analysis.atmosphere.values.end()) {
          const real mismatch =
              std::fabs(geopotential_height->second[idx] -
                        metrics.z_center(i, j, k));
          if (diagnostics.max_abs_z_src_minus_metric.has_value()) {
            diagnostics.max_abs_z_src_minus_metric = std::max(
                *diagnostics.max_abs_z_src_minus_metric, mismatch);
          } else {
            diagnostics.max_abs_z_src_minus_metric = mismatch;
          }
        }
      }

      for (int k_face = 1; k_face < analysis.grid.nz; ++k_face) {
        const int k_lower = k_face - 1;
        const int k_upper = k_face;
        const auto idx_lower =
            linear_index_3d(i, j, k_lower, analysis.grid.nx, analysis.grid.ny);
        const auto idx_upper =
            linear_index_3d(i, j, k_upper, analysis.grid.nx, analysis.grid.ny);
        const real p_pert_lower =
            pressure_global[idx_lower] - balanced.pressure[idx_lower];
        const real p_pert_upper =
            pressure_global[idx_upper] - balanced.pressure[idx_upper];
        const real rho_pert_face =
            0.5f * ((rho_global[idx_lower] - balanced.rho_d[idx_lower]) +
                    (rho_global[idx_upper] - balanced.rho_d[idx_upper]));
        const real rho_ref_face =
            0.5f * (balanced.rho_d[idx_lower] + balanced.rho_d[idx_upper]);
        const real scale = std::max(
            {std::fabs(p_pert_lower) /
                 std::max(std::fabs(balanced.pressure[idx_lower]), 1.0f),
             std::fabs(p_pert_upper) /
                 std::max(std::fabs(balanced.pressure[idx_upper]), 1.0f),
             std::fabs(rho_pert_face) / std::max(std::fabs(rho_ref_face), 1.0f),
             0.0f});
        diagnostics.max_rel_fast_vertical = std::max(
            diagnostics.max_rel_fast_vertical, scale);
      }
    }
  }

  if (tracers != nullptr) {
    for (std::size_t rank = 0; rank < layout.size(); ++rank) {
      const auto& desc = layout[rank];
      const auto& tracer_state = tracers->at(rank);
      for (std::size_t tracer_index = 0; tracer_index < tracer_state.size();
           ++tracer_index) {
        const auto& spec = tracer_state.registry().at(static_cast<int>(tracer_index));
        const auto field_it = analysis.atmosphere.values.find(spec.name);
        if (field_it == analysis.atmosphere.values.end()) {
          continue;
        }

        const auto& mixing_ratio = field_it->second;
        const auto& rho_q = tracer_state.mass(static_cast<int>(tracer_index));
        for (int j = 0; j < desc.ny_local(); ++j) {
          const int j_global = desc.j_begin + j;
          for (int i = 0; i < desc.nx_local(); ++i) {
            const int i_global = desc.i_begin + i;
            for (int k = 0; k < desc.nz; ++k) {
              const auto idx = linear_index_3d(i_global, j_global, k,
                                               analysis.grid.nx,
                                               analysis.grid.ny);
              const real expected =
                  dry_states[rank].rho_d(i, j, k) * mixing_ratio[idx];
              diagnostics.max_abs_tracer_closure = std::max(
                  diagnostics.max_abs_tracer_closure,
                  std::fabs(rho_q(i, j, k) - expected));
            }
          }
        }
      }
    }
  }

  return diagnostics;
}

std::vector<state::TracerState> make_specific_humidity_tracers_from_analysis(
    const AnalysisStateIR& analysis,
    const std::vector<dycore::DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const std::string& label_prefix) {
  return make_tracers_from_analysis(
      analysis, dry_states, layout, state::make_specific_humidity_registry(),
      label_prefix);
}

std::vector<state::TracerState> make_tracers_from_analysis(
    const AnalysisStateIR& analysis,
    const std::vector<dycore::DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const state::TracerRegistry& registry, const std::string& label_prefix,
    const std::vector<state::TracerState>* baseline_tracers) {
  validate_runtime_analysis(analysis);
  gwm::require(dry_states.size() == layout.size(),
               "Dry-state/layout size mismatch in make_tracers_from_analysis");
  if (baseline_tracers != nullptr) {
    gwm::require(baseline_tracers->size() == layout.size(),
                 "Baseline tracer/layout size mismatch in "
                 "make_tracers_from_analysis");
  }

  std::vector<state::TracerState> tracers;
  tracers.reserve(layout.size());
  for (std::size_t rank = 0; rank < layout.size(); ++rank) {
    const auto& desc = layout[rank];
    state::TracerState tracer_state(registry, desc.nx_local(), desc.ny_local(),
                                    desc.nz, desc.halo,
                                    label_prefix + "_rank_" +
                                        std::to_string(desc.rank));
    tracer_state.fill_zero();
    if (baseline_tracers != nullptr) {
      tracer_state.copy_all_from(baseline_tracers->at(rank));
    }
    tracers.push_back(std::move(tracer_state));
  }

  for (std::size_t tracer_index = 0; tracer_index < registry.size();
       ++tracer_index) {
    const auto& spec = registry.at(static_cast<int>(tracer_index));
    if (spec.name == state::kSpecificHumidityTracerName) {
      const auto& specific_humidity =
          analysis.atmosphere.values.at(state::kSpecificHumidityTracerName);
      gwm::require(
          specific_humidity.size() == analysis.grid.cell_count_3d(),
          "specific_humidity analysis field size mismatch");
      for (std::size_t rank = 0; rank < layout.size(); ++rank) {
        const auto& desc = layout[rank];
        auto& rho_q = tracers[rank].mass(spec.name);
        for (int j = 0; j < desc.ny_local(); ++j) {
          const int j_global = desc.j_begin + j;
          for (int i = 0; i < desc.nx_local(); ++i) {
            const int i_global = desc.i_begin + i;
            for (int k = 0; k < desc.nz; ++k) {
              const auto idx = linear_index_3d(i_global, j_global, k,
                                               analysis.grid.nx,
                                               analysis.grid.ny);
              rho_q(i, j, k) =
                  dry_states[rank].rho_d(i, j, k) * specific_humidity[idx];
            }
          }
        }
      }
      continue;
    }

    const auto field_it = analysis.atmosphere.values.find(spec.name);
    if (field_it == analysis.atmosphere.values.end()) {
      continue;
    }
    const auto& mixing_ratio = field_it->second;
    gwm::require(
        mixing_ratio.size() == analysis.grid.cell_count_3d(),
        "Tracer analysis field size mismatch for " + spec.name);
    for (std::size_t rank = 0; rank < layout.size(); ++rank) {
      const auto& desc = layout[rank];
      auto& rho_q = tracers[rank].mass(spec.name);
      for (int j = 0; j < desc.ny_local(); ++j) {
        const int j_global = desc.j_begin + j;
        for (int i = 0; i < desc.nx_local(); ++i) {
          const int i_global = desc.i_begin + i;
          for (int k = 0; k < desc.nz; ++k) {
            const auto idx = linear_index_3d(i_global, j_global, k,
                                             analysis.grid.nx,
                                             analysis.grid.ny);
            rho_q(i, j, k) =
                dry_states[rank].rho_d(i, j, k) * mixing_ratio[idx];
          }
        }
      }
    }
  }

  return tracers;
}

std::vector<state::TracerState> make_warm_rain_tracers_from_analysis(
    const AnalysisStateIR& analysis,
    const std::vector<dycore::DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const std::string& label_prefix) {
  return make_tracers_from_analysis(analysis, dry_states, layout,
                                    state::make_warm_rain_registry(),
                                    label_prefix);
}

PreparedCaseBoundaryUpdater::PreparedCaseBoundaryUpdater(
    AnalysisStateIR analysis, BoundaryCacheIR cache, PreparedCaseInitConfig config)
    : analysis_(std::move(analysis)),
      cache_(std::move(cache)),
      config_(config) {
  metrics_ = make_prepared_case_metrics(analysis_, config_);
}

void PreparedCaseBoundaryUpdater::set_step_start_time(
    real step_start_time_seconds) {
  step_start_time_seconds_ = std::max(step_start_time_seconds, 0.0f);
}

void PreparedCaseBoundaryUpdater::apply(
    std::vector<dycore::DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout, real sim_time) {
  if (cache_.empty()) {
    return;
  }

  refresh_boundary_projection_cache(
      analysis_, cache_, metrics_, step_start_time_seconds_, sim_time,
      cached_forecast_offset_seconds_, cached_snapshot_, cached_rho_d_,
      cached_rho_theta_m_, cached_mom_u_, cached_mom_v_, cached_mom_w_);
  apply_cached_dry_boundary_strips(states, layout, analysis_.grid.nx,
                                   analysis_.grid.ny, cached_rho_d_,
                                   cached_rho_theta_m_, cached_mom_u_,
                                   cached_mom_v_, cached_mom_w_);
}

PreparedCaseTracerBoundaryUpdater::PreparedCaseTracerBoundaryUpdater(
    AnalysisStateIR analysis, BoundaryCacheIR cache, PreparedCaseInitConfig config)
    : analysis_(std::move(analysis)),
      cache_(std::move(cache)),
      config_(config) {
  metrics_ = make_prepared_case_metrics(analysis_, config_);
}

void PreparedCaseTracerBoundaryUpdater::set_step_start_time(
    real step_start_time_seconds) {
  step_start_time_seconds_ = std::max(step_start_time_seconds, 0.0f);
}

void PreparedCaseTracerBoundaryUpdater::apply(
    std::vector<state::TracerState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout, real sim_time) {
  if (cache_.empty() || states.empty()) {
    return;
  }

  refresh_boundary_projection_cache(
      analysis_, cache_, metrics_, step_start_time_seconds_, sim_time,
      cached_forecast_offset_seconds_, cached_snapshot_, cached_rho_d_,
      cached_rho_theta_m_, cached_mom_u_, cached_mom_v_, cached_mom_w_);
  apply_cached_tracer_boundary_strips(states, layout, cached_snapshot_,
                                      analysis_.grid.nx, analysis_.grid.ny,
                                      cached_rho_d_);
}

}  // namespace gwm::ingest
