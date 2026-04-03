#include "gwm/dycore/dry_fast_modes.hpp"

#include <algorithm>
#include <vector>

#include "gwm/comm/cartesian_topology.hpp"
#include "gwm/comm/halo_exchange.hpp"
#include "gwm/core/cuda_utils.hpp"
#include "gwm/core/dry_thermo.hpp"
#include "gwm/dycore/dry_pressure_gradient.hpp"

namespace gwm::dycore {

namespace {

__host__ __device__ std::size_t metric_storage_index(int i, int j, int k,
                                                     int nx, int ny) {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

__host__ __device__ int wrap_metric_index(int idx, int n, bool periodic) {
  if (!periodic) {
    return idx;
  }
  const int mod = idx % n;
  return mod < 0 ? mod + n : mod;
}

__host__ __device__ std::size_t face_storage_index(int i, int j, int k, int nx,
                                                   int ny, int halo) {
  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  return (static_cast<std::size_t>(k) * pitch_y +
          static_cast<std::size_t>(j + halo)) *
             static_cast<std::size_t>(pitch_x) +
         static_cast<std::size_t>(i + halo);
}

__global__ void accumulate_level_reference_kernel(
    const real* rho_d, const real* pressure, real* rho_sum, real* pressure_sum,
    int nx, int ny, int nz, int halo) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  const int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (i >= nx || j >= ny || k >= nz) {
    return;
  }

  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  const auto cell_idx = (k * pitch_y + (j + halo)) * pitch_x + (i + halo);
  atomicAdd(&rho_sum[k], rho_d[cell_idx]);
  atomicAdd(&pressure_sum[k], pressure[cell_idx]);
}

__global__ void finalize_level_reference_kernel(real* rho_ref, real* pressure_ref,
                                                int nz, real inv_count) {
  const int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= nz) {
    return;
  }
  rho_ref[k] *= inv_count;
  pressure_ref[k] *= inv_count;
}

__global__ void fill_reference_from_levels_kernel(
    real* rho_ref_field, real* pressure_ref_field, const real* rho_ref_levels,
    const real* pressure_ref_levels, int nx, int ny, int nz, int halo) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  const int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (i >= nx || j >= ny || k >= nz) {
    return;
  }

  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  const auto cell_idx = (k * pitch_y + (j + halo)) * pitch_x + (i + halo);
  rho_ref_field[cell_idx] = rho_ref_levels[k];
  pressure_ref_field[cell_idx] = pressure_ref_levels[k];
}

__global__ void build_column_hydrostatic_reference_kernel(
    const real* rho_d, const real* rho_theta_m, const real* pressure,
    const real* z_centers_metric, real* rho_ref, real* pressure_ref,
    int metrics_nx, int metrics_ny, bool periodic_x, bool periodic_y,
    int i_begin, int j_begin, int nx, int ny, int nz, int halo, real gravity) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  if (i >= nx || j >= ny) {
    return;
  }

  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  const auto cell_idx = [&](int ii, int jj, int kk) {
    return (kk * pitch_y + (jj + halo)) * pitch_x + (ii + halo);
  };
  const auto metric_idx = [&](int kk) {
    return metric_storage_index(
        wrap_metric_index(i_begin + i, metrics_nx, periodic_x),
        wrap_metric_index(j_begin + j, metrics_ny, periodic_y), kk, metrics_nx,
        metrics_ny);
  };

  const real theta0 = gwm::core::dry_theta_from_conserved(
      rho_d[cell_idx(i, j, 0)], rho_theta_m[cell_idx(i, j, 0)]);
  pressure_ref[cell_idx(i, j, 0)] = pressure[cell_idx(i, j, 0)];
  rho_ref[cell_idx(i, j, 0)] = gwm::core::dry_rho_from_pressure_theta(
      pressure_ref[cell_idx(i, j, 0)], theta0);

  for (int k = 1; k < nz; ++k) {
    const real theta_k = gwm::core::dry_theta_from_conserved(
        rho_d[cell_idx(i, j, k)], rho_theta_m[cell_idx(i, j, k)]);
    const real dz =
        z_centers_metric[metric_idx(k)] - z_centers_metric[metric_idx(k - 1)];

    real p_next = fmaxf(1.0f, pressure_ref[cell_idx(i, j, k - 1)] -
                                  gravity * rho_ref[cell_idx(i, j, k - 1)] * dz);
    for (int iter = 0; iter < 4; ++iter) {
      const real rho_next =
          gwm::core::dry_rho_from_pressure_theta(p_next, theta_k);
      p_next = pressure_ref[cell_idx(i, j, k - 1)] -
               gravity * 0.5f * (rho_ref[cell_idx(i, j, k - 1)] + rho_next) *
                   dz;
      p_next = fmaxf(p_next, 1.0f);
    }

    pressure_ref[cell_idx(i, j, k)] = p_next;
    rho_ref[cell_idx(i, j, k)] =
        gwm::core::dry_rho_from_pressure_theta(p_next, theta_k);
  }
}

__global__ void update_fast_x_momentum_kernel(
    const real* pressure, const real* z_centers_metric, real* mom_u,
    int metrics_nx, int metrics_ny, bool periodic_x, bool periodic_y,
    int i_begin, int j_begin, int nx, int ny, int nz, int halo, real dx,
    real dt_fast, bool open_west, bool open_east) {
  const int i_face = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  const int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (i_face > nx || j >= ny || k >= nz) {
    return;
  }

  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  const auto cell_idx = [&](int ii, int jj, int kk) {
    return (kk * pitch_y + (jj + halo)) * pitch_x + (ii + halo);
  };
  const auto face_idx = [&](int ii, int jj, int kk) {
    const int face_pitch_x = (nx + 1) + 2 * halo;
    const int face_pitch_y = ny + 2 * halo;
    return (kk * face_pitch_y + (jj + halo)) * face_pitch_x + (ii + halo);
  };

  if ((i_face == 0 && open_west) || (i_face == nx && open_east)) {
    return;
  }

  const real p_left = pressure[cell_idx(i_face - 1, j, k)];
  const real p_right = pressure[cell_idx(i_face, j, k)];
  const int gj = wrap_metric_index(j_begin + j, metrics_ny, periodic_y);
  const int gi_left =
      wrap_metric_index(i_begin + i_face - 1, metrics_nx, periodic_x);
  const int gi_right =
      wrap_metric_index(i_begin + i_face, metrics_nx, periodic_x);
  const auto vertical_gradient = [&](int gi, int local_i) {
    if (nz == 1) {
      return 0.0f;
    }
    if (k == 0) {
      const real dz =
          z_centers_metric[metric_storage_index(gi, gj, 1, metrics_nx, metrics_ny)] -
          z_centers_metric[metric_storage_index(gi, gj, 0, metrics_nx, metrics_ny)];
      return (pressure[cell_idx(local_i, j, 1)] - pressure[cell_idx(local_i, j, 0)]) /
             fmaxf(dz, 1.0e-6f);
    }
    if (k + 1 == nz) {
      const real dz =
          z_centers_metric[metric_storage_index(gi, gj, k, metrics_nx, metrics_ny)] -
          z_centers_metric[metric_storage_index(gi, gj, k - 1, metrics_nx, metrics_ny)];
      return (pressure[cell_idx(local_i, j, k)] -
              pressure[cell_idx(local_i, j, k - 1)]) /
             fmaxf(dz, 1.0e-6f);
    }
    const real dz =
        z_centers_metric[metric_storage_index(gi, gj, k + 1, metrics_nx, metrics_ny)] -
        z_centers_metric[metric_storage_index(gi, gj, k - 1, metrics_nx, metrics_ny)];
    return (pressure[cell_idx(local_i, j, k + 1)] -
            pressure[cell_idx(local_i, j, k - 1)]) /
           fmaxf(dz, 1.0e-6f);
  };
  const real z_left =
      z_centers_metric[metric_storage_index(gi_left, gj, k, metrics_nx, metrics_ny)];
  const real z_right =
      z_centers_metric[metric_storage_index(gi_right, gj, k, metrics_nx, metrics_ny)];
  const real dzdx = (z_right - z_left) / dx;
  const real dpdz_face = 0.5f * (vertical_gradient(gi_left, i_face - 1) +
                                 vertical_gradient(gi_right, i_face));
  mom_u[face_idx(i_face, j, k)] +=
      -dt_fast * ((p_right - p_left) / dx - dzdx * dpdz_face);
}

__global__ void update_fast_y_momentum_kernel(
    const real* pressure, const real* z_centers_metric, real* mom_v,
    int metrics_nx, int metrics_ny, bool periodic_x, bool periodic_y,
    int i_begin, int j_begin, int nx, int ny, int nz, int halo, real dy,
    real dt_fast, bool open_south, bool open_north) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j_face = blockIdx.y * blockDim.y + threadIdx.y;
  const int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (i >= nx || j_face > ny || k >= nz) {
    return;
  }

  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  const auto cell_idx = [&](int ii, int jj, int kk) {
    return (kk * pitch_y + (jj + halo)) * pitch_x + (ii + halo);
  };
  const auto face_idx = [&](int ii, int jj, int kk) {
    const int face_pitch_x = nx + 2 * halo;
    const int face_pitch_y = (ny + 1) + 2 * halo;
    return (kk * face_pitch_y + (jj + halo)) * face_pitch_x + (ii + halo);
  };

  if ((j_face == 0 && open_south) || (j_face == ny && open_north)) {
    return;
  }

  const real p_south = pressure[cell_idx(i, j_face - 1, k)];
  const real p_north = pressure[cell_idx(i, j_face, k)];
  const int gi = wrap_metric_index(i_begin + i, metrics_nx, periodic_x);
  const int gj_south =
      wrap_metric_index(j_begin + j_face - 1, metrics_ny, periodic_y);
  const int gj_north =
      wrap_metric_index(j_begin + j_face, metrics_ny, periodic_y);
  const auto vertical_gradient = [&](int gj, int local_j) {
    if (nz == 1) {
      return 0.0f;
    }
    if (k == 0) {
      const real dz =
          z_centers_metric[metric_storage_index(gi, gj, 1, metrics_nx, metrics_ny)] -
          z_centers_metric[metric_storage_index(gi, gj, 0, metrics_nx, metrics_ny)];
      return (pressure[cell_idx(i, local_j, 1)] - pressure[cell_idx(i, local_j, 0)]) /
             fmaxf(dz, 1.0e-6f);
    }
    if (k + 1 == nz) {
      const real dz =
          z_centers_metric[metric_storage_index(gi, gj, k, metrics_nx, metrics_ny)] -
          z_centers_metric[metric_storage_index(gi, gj, k - 1, metrics_nx, metrics_ny)];
      return (pressure[cell_idx(i, local_j, k)] -
              pressure[cell_idx(i, local_j, k - 1)]) /
             fmaxf(dz, 1.0e-6f);
    }
    const real dz =
        z_centers_metric[metric_storage_index(gi, gj, k + 1, metrics_nx, metrics_ny)] -
        z_centers_metric[metric_storage_index(gi, gj, k - 1, metrics_nx, metrics_ny)];
    return (pressure[cell_idx(i, local_j, k + 1)] -
            pressure[cell_idx(i, local_j, k - 1)]) /
           fmaxf(dz, 1.0e-6f);
  };
  const real z_south =
      z_centers_metric[metric_storage_index(gi, gj_south, k, metrics_nx, metrics_ny)];
  const real z_north =
      z_centers_metric[metric_storage_index(gi, gj_north, k, metrics_nx, metrics_ny)];
  const real dzdy = (z_north - z_south) / dy;
  const real dpdz_face =
      0.5f * (vertical_gradient(gj_south, j_face - 1) +
              vertical_gradient(gj_north, j_face));
  mom_v[face_idx(i, j_face, k)] +=
      -dt_fast * ((p_north - p_south) / dy - dzdy * dpdz_face);
}

__global__ void update_fast_vertical_momentum_kernel(
    const real* pressure, const real* rho_d, const real* pressure_ref,
    const real* rho_ref, const real* inv_dz_face_metric, real* mom_w,
    int metrics_nx, int metrics_ny, bool periodic_x, bool periodic_y,
    int i_begin, int j_begin, int nx, int ny, int nz, int halo, real dt_fast,
    real gravity) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  const int k_face = blockIdx.z * blockDim.z + threadIdx.z;
  if (i >= nx || j >= ny || k_face > nz) {
    return;
  }

  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  const auto cell_idx = [&](int ii, int jj, int kk) {
    return (kk * pitch_y + (jj + halo)) * pitch_x + (ii + halo);
  };
  const auto zface_idx = [&](int ii, int jj, int kk) {
    return face_storage_index(ii, jj, kk, nx, ny, halo);
  };

  if (k_face == 0 || k_face == nz) {
    return;
  }

  const int k_lower = k_face - 1;
  const int k_upper = k_face;
  const real p_pert_lower =
      pressure[cell_idx(i, j, k_lower)] - pressure_ref[cell_idx(i, j, k_lower)];
  const real p_pert_upper =
      pressure[cell_idx(i, j, k_upper)] - pressure_ref[cell_idx(i, j, k_upper)];
  const real rho_pert_face =
      0.5f * ((rho_d[cell_idx(i, j, k_lower)] -
               rho_ref[cell_idx(i, j, k_lower)]) +
              (rho_d[cell_idx(i, j, k_upper)] -
               rho_ref[cell_idx(i, j, k_upper)]));
  const real inv_dz_face =
      inv_dz_face_metric[metric_storage_index(
          wrap_metric_index(i_begin + i, metrics_nx, periodic_x),
          wrap_metric_index(j_begin + j, metrics_ny, periodic_y), k_face,
          metrics_nx, metrics_ny)];
  const real tendency =
      -(p_pert_upper - p_pert_lower) * inv_dz_face - gravity * rho_pert_face;
  mom_w[zface_idx(i, j, k_face)] += dt_fast * tendency;
}

__global__ void update_density_from_divergence_kernel(
    real* rho_d, const real* mom_u, const real* mom_v, const real* mom_w,
    const real* inv_dz_cell_metric, int metrics_nx, int metrics_ny,
    bool periodic_x, bool periodic_y, int i_begin, int j_begin, int nx, int ny,
    int nz, int halo, real dx, real dy, real dt_fast) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  const int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (i >= nx || j >= ny || k >= nz) {
    return;
  }

  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  const auto cell_idx = [&](int ii, int jj, int kk) {
    return (kk * pitch_y + (jj + halo)) * pitch_x + (ii + halo);
  };
  const auto xface_idx = [&](int ii, int jj, int kk) {
    return face_storage_index(ii, jj, kk, nx + 1, ny, halo);
  };
  const auto yface_idx = [&](int ii, int jj, int kk) {
    return face_storage_index(ii, jj, kk, nx, ny + 1, halo);
  };
  const auto zface_idx = [&](int ii, int jj, int kk) {
    return face_storage_index(ii, jj, kk, nx, ny, halo);
  };

  const real div_mass =
      (mom_u[xface_idx(i + 1, j, k)] - mom_u[xface_idx(i, j, k)]) / dx +
      (mom_v[yface_idx(i, j + 1, k)] - mom_v[yface_idx(i, j, k)]) / dy +
      (mom_w[zface_idx(i, j, k + 1)] - mom_w[zface_idx(i, j, k)]) *
          inv_dz_cell_metric[metric_storage_index(
              wrap_metric_index(i_begin + i, metrics_nx, periodic_x),
              wrap_metric_index(j_begin + j, metrics_ny, periodic_y), k,
              metrics_nx, metrics_ny)];

  rho_d[cell_idx(i, j, k)] =
      fmaxf(1.0e-6f, rho_d[cell_idx(i, j, k)] - dt_fast * div_mass);
}

void sync_after_launches() { GWM_CUDA_CHECK(cudaDeviceSynchronize()); }

}  // namespace

void apply_local_split_explicit_fast_modes(
    std::vector<DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const DryStepperConfig& config) {
  gwm::require(states.size() == layout.size(),
               "State/layout size mismatch in apply_local_split_explicit_fast_modes");
  if (states.empty()) {
    return;
  }

  const int substeps = std::max(1, config.fast_substeps);
  const real dt_fast = config.dt / static_cast<real>(substeps);
  std::vector<state::Field3D<real>> rho_theta_halo;
  std::vector<state::Field3D<real>> rho_halo;
  std::vector<state::Field3D<real>> pressure_fields;
  std::vector<state::Field3D<real>> rho_ref_fields;
  std::vector<state::Field3D<real>> pressure_ref_fields;
  rho_theta_halo.reserve(states.size());
  rho_halo.reserve(states.size());
  pressure_fields.reserve(states.size());
  rho_ref_fields.reserve(states.size());
  pressure_ref_fields.reserve(states.size());
  for (const auto& state : states) {
    rho_theta_halo.push_back(state.rho_theta_m.clone_empty_like("_fast_theta"));
    rho_theta_halo.back().copy_all_from(state.rho_theta_m);
    rho_halo.push_back(state.rho_d.clone_empty_like("_fast_rho"));
    pressure_fields.push_back(state.rho_d.clone_empty_like("_fast_pressure"));
    rho_ref_fields.push_back(state.rho_d.clone_empty_like("_fast_rho_ref"));
    pressure_ref_fields.push_back(
        state.rho_d.clone_empty_like("_fast_pressure_ref"));
  }
  comm::HaloExchange::exchange_scalar(rho_theta_halo, layout);

  real* rho_ref_levels = nullptr;
  real* pressure_ref_levels = nullptr;
  GWM_CUDA_CHECK(cudaMallocManaged(&rho_ref_levels, sizeof(real) * metrics.nz));
  GWM_CUDA_CHECK(
      cudaMallocManaged(&pressure_ref_levels, sizeof(real) * metrics.nz));

  int horizontal_cell_count = 0;
  for (const auto& desc : layout) {
    horizontal_cell_count += desc.nx_local() * desc.ny_local();
  }
  const real inv_horizontal_count =
      1.0f / static_cast<real>(std::max(1, horizontal_cell_count));

  bool terrain_active = false;
  for (int j = 0; j < metrics.ny && !terrain_active; ++j) {
    for (int i = 0; i < metrics.nx; ++i) {
      if (std::abs(static_cast<double>(metrics.terrain(i, j, 0))) > 1.0e-6) {
        terrain_active = true;
        break;
      }
    }
  }

  dim3 block(8, 8, 4);
  for (int substep = 0; substep < substeps; ++substep) {
    for (std::size_t n = 0; n < states.size(); ++n) {
      rho_halo[n].copy_all_from(states[n].rho_d);
    }
    comm::HaloExchange::exchange_scalar(rho_halo, layout);
    compute_dry_pressure_fields(rho_halo, rho_theta_halo, pressure_fields);

    if (!terrain_active) {
      GWM_CUDA_CHECK(cudaMemset(rho_ref_levels, 0, sizeof(real) * metrics.nz));
      GWM_CUDA_CHECK(
          cudaMemset(pressure_ref_levels, 0, sizeof(real) * metrics.nz));

      for (std::size_t n = 0; n < states.size(); ++n) {
        dim3 cell_grid((states[n].rho_d.nx() + block.x - 1) / block.x,
                       (states[n].rho_d.ny() + block.y - 1) / block.y,
                       (states[n].rho_d.nz() + block.z - 1) / block.z);
        accumulate_level_reference_kernel<<<cell_grid, block>>>(
            rho_halo[n].data(), pressure_fields[n].data(), rho_ref_levels,
            pressure_ref_levels, states[n].rho_d.nx(), states[n].rho_d.ny(),
            states[n].rho_d.nz(), states[n].rho_d.halo());
        GWM_CUDA_CHECK(cudaGetLastError());
      }

      constexpr int level_block = 256;
      const int level_grid = (metrics.nz + level_block - 1) / level_block;
      finalize_level_reference_kernel<<<level_grid, level_block>>>(
          rho_ref_levels, pressure_ref_levels, metrics.nz, inv_horizontal_count);
      GWM_CUDA_CHECK(cudaGetLastError());
      sync_after_launches();
    }

    std::vector<state::FaceField<real>*> u_faces;
    std::vector<state::FaceField<real>*> v_faces;
    std::vector<state::FaceField<real>*> w_faces;
    u_faces.reserve(states.size());
    v_faces.reserve(states.size());
    w_faces.reserve(states.size());

    for (std::size_t n = 0; n < states.size(); ++n) {
      u_faces.push_back(&states[n].mom_u);
      v_faces.push_back(&states[n].mom_v);
      w_faces.push_back(&states[n].mom_w);

      const auto neighbors = comm::find_cartesian_neighbors(layout, layout[n]);
      if (terrain_active) {
        dim3 column_grid((states[n].rho_d.nx() + block.x - 1) / block.x,
                         (states[n].rho_d.ny() + block.y - 1) / block.y, 1);
        build_column_hydrostatic_reference_kernel<<<column_grid, dim3(block.x, block.y, 1)>>>(
            rho_halo[n].data(), rho_theta_halo[n].data(),
            pressure_fields[n].data(), metrics.z_centers_metric.data(),
            rho_ref_fields[n].data(), pressure_ref_fields[n].data(),
            metrics.nx, metrics.ny, metrics.periodic_x, metrics.periodic_y,
            layout[n].i_begin, layout[n].j_begin, states[n].rho_d.nx(),
            states[n].rho_d.ny(), states[n].rho_d.nz(),
            states[n].rho_d.halo(), config.gravity);
      } else {
        dim3 cell_grid((states[n].rho_d.nx() + block.x - 1) / block.x,
                       (states[n].rho_d.ny() + block.y - 1) / block.y,
                       (states[n].rho_d.nz() + block.z - 1) / block.z);
        fill_reference_from_levels_kernel<<<cell_grid, block>>>(
            rho_ref_fields[n].data(), pressure_ref_fields[n].data(),
            rho_ref_levels, pressure_ref_levels, states[n].rho_d.nx(),
            states[n].rho_d.ny(), states[n].rho_d.nz(),
            states[n].rho_d.halo());
      }
      GWM_CUDA_CHECK(cudaGetLastError());

      dim3 xgrid((states[n].rho_d.nx() + 1 + block.x - 1) / block.x,
                 (states[n].rho_d.ny() + block.y - 1) / block.y,
                 (states[n].rho_d.nz() + block.z - 1) / block.z);
      update_fast_x_momentum_kernel<<<xgrid, block>>>(
          pressure_fields[n].data(), metrics.z_centers_metric.data(),
          states[n].mom_u.storage().data(), metrics.nx, metrics.ny,
          metrics.periodic_x, metrics.periodic_y,
          layout[n].i_begin, layout[n].j_begin, states[n].rho_d.nx(),
          states[n].rho_d.ny(), states[n].rho_d.nz(), states[n].rho_d.halo(),
          static_cast<real>(metrics.dx), dt_fast, neighbors.west < 0,
          neighbors.east < 0);
      GWM_CUDA_CHECK(cudaGetLastError());

      dim3 ygrid((states[n].rho_d.nx() + block.x - 1) / block.x,
                 (states[n].rho_d.ny() + 1 + block.y - 1) / block.y,
                 (states[n].rho_d.nz() + block.z - 1) / block.z);
      update_fast_y_momentum_kernel<<<ygrid, block>>>(
          pressure_fields[n].data(), metrics.z_centers_metric.data(),
          states[n].mom_v.storage().data(), metrics.nx, metrics.ny,
          metrics.periodic_x, metrics.periodic_y,
          layout[n].i_begin, layout[n].j_begin, states[n].rho_d.nx(),
          states[n].rho_d.ny(), states[n].rho_d.nz(), states[n].rho_d.halo(),
          static_cast<real>(metrics.dy), dt_fast, neighbors.south < 0,
          neighbors.north < 0);
      GWM_CUDA_CHECK(cudaGetLastError());

      dim3 zgrid((states[n].rho_d.nx() + block.x - 1) / block.x,
                 (states[n].rho_d.ny() + block.y - 1) / block.y,
                 (states[n].rho_d.nz() + 1 + block.z - 1) / block.z);
      update_fast_vertical_momentum_kernel<<<zgrid, block>>>(
          pressure_fields[n].data(), rho_halo[n].data(),
          pressure_ref_fields[n].data(), rho_ref_fields[n].data(),
          metrics.inv_dz_face_metric.data(), states[n].mom_w.storage().data(),
          metrics.nx, metrics.ny, metrics.periodic_x, metrics.periodic_y,
          layout[n].i_begin, layout[n].j_begin, states[n].rho_d.nx(),
          states[n].rho_d.ny(), states[n].rho_d.nz(), states[n].rho_d.halo(),
          dt_fast, config.gravity);
      GWM_CUDA_CHECK(cudaGetLastError());
    }
    sync_after_launches();

    comm::HaloExchange::exchange_face(u_faces, layout);
    comm::HaloExchange::exchange_face(v_faces, layout);
    comm::HaloExchange::exchange_face(w_faces, layout);
    comm::HaloExchange::synchronize_owned_face_interfaces(u_faces, layout);
    comm::HaloExchange::synchronize_owned_face_interfaces(v_faces, layout);
    comm::HaloExchange::exchange_face(u_faces, layout);
    comm::HaloExchange::exchange_face(v_faces, layout);

    for (std::size_t n = 0; n < states.size(); ++n) {
      dim3 cell_grid((states[n].rho_d.nx() + block.x - 1) / block.x,
                     (states[n].rho_d.ny() + block.y - 1) / block.y,
                     (states[n].rho_d.nz() + block.z - 1) / block.z);
      update_density_from_divergence_kernel<<<cell_grid, block>>>(
          states[n].rho_d.data(), states[n].mom_u.storage().data(),
          states[n].mom_v.storage().data(), states[n].mom_w.storage().data(),
          metrics.inv_dz_cell_metric.data(), metrics.nx, metrics.ny,
          metrics.periodic_x, metrics.periodic_y, layout[n].i_begin,
          layout[n].j_begin, states[n].rho_d.nx(), states[n].rho_d.ny(),
          states[n].rho_d.nz(), states[n].rho_d.halo(),
          static_cast<real>(metrics.dx), static_cast<real>(metrics.dy),
          dt_fast);
      GWM_CUDA_CHECK(cudaGetLastError());
    }
    sync_after_launches();
  }

  GWM_CUDA_CHECK(cudaFree(rho_ref_levels));
  GWM_CUDA_CHECK(cudaFree(pressure_ref_levels));
}

}  // namespace gwm::dycore
