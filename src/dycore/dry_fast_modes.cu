#include "gwm/dycore/dry_fast_modes.hpp"

#include <algorithm>
#include <vector>

#include "gwm/comm/cartesian_topology.hpp"
#include "gwm/comm/halo_exchange.hpp"
#include "gwm/core/cuda_utils.hpp"
#include "gwm/dycore/dry_pressure_gradient.hpp"

namespace gwm::dycore {

namespace {

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

__global__ void update_fast_x_momentum_kernel(
    const real* pressure, real* mom_u, int nx, int ny, int nz, int halo,
    real dx, real dt_fast, bool open_west, bool open_east) {
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
  mom_u[face_idx(i_face, j, k)] += -dt_fast * (p_right - p_left) / dx;
}

__global__ void update_fast_y_momentum_kernel(
    const real* pressure, real* mom_v, int nx, int ny, int nz, int halo,
    real dy, real dt_fast, bool open_south, bool open_north) {
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
  mom_v[face_idx(i, j_face, k)] += -dt_fast * (p_north - p_south) / dy;
}

__global__ void update_fast_vertical_momentum_kernel(
    const real* pressure, const real* rho_d, const real* pressure_ref,
    const real* rho_ref, real* mom_w, int nx, int ny, int nz, int halo,
    real dz, real dt_fast, real gravity) {
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
      pressure[cell_idx(i, j, k_lower)] - pressure_ref[k_lower];
  const real p_pert_upper =
      pressure[cell_idx(i, j, k_upper)] - pressure_ref[k_upper];
  const real rho_pert_face =
      0.5f * ((rho_d[cell_idx(i, j, k_lower)] - rho_ref[k_lower]) +
              (rho_d[cell_idx(i, j, k_upper)] - rho_ref[k_upper]));
  const real tendency =
      -(p_pert_upper - p_pert_lower) / dz - gravity * rho_pert_face;
  mom_w[zface_idx(i, j, k_face)] += dt_fast * tendency;
}

__global__ void update_density_from_divergence_kernel(
    real* rho_d, const real* mom_u, const real* mom_v, const real* mom_w,
    int nx, int ny, int nz, int halo, real dx, real dy, real dz, real dt_fast) {
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
      (mom_w[zface_idx(i, j, k + 1)] - mom_w[zface_idx(i, j, k)]) / dz;

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
  const real dz = static_cast<real>(metrics.dz_nominal);

  std::vector<state::Field3D<real>> rho_theta_halo;
  std::vector<state::Field3D<real>> rho_halo;
  std::vector<state::Field3D<real>> pressure_fields;
  rho_theta_halo.reserve(states.size());
  rho_halo.reserve(states.size());
  pressure_fields.reserve(states.size());
  for (const auto& state : states) {
    rho_theta_halo.push_back(state.rho_theta_m.clone_empty_like("_fast_theta"));
    rho_theta_halo.back().copy_all_from(state.rho_theta_m);
    rho_halo.push_back(state.rho_d.clone_empty_like("_fast_rho"));
    pressure_fields.push_back(state.rho_d.clone_empty_like("_fast_pressure"));
  }
  comm::HaloExchange::exchange_scalar(rho_theta_halo, layout);

  real* rho_ref = nullptr;
  real* pressure_ref = nullptr;
  GWM_CUDA_CHECK(cudaMallocManaged(&rho_ref, sizeof(real) * metrics.nz));
  GWM_CUDA_CHECK(cudaMallocManaged(&pressure_ref, sizeof(real) * metrics.nz));

  int horizontal_cell_count = 0;
  for (const auto& desc : layout) {
    horizontal_cell_count += desc.nx_local() * desc.ny_local();
  }
  const real inv_horizontal_count =
      1.0f / static_cast<real>(std::max(1, horizontal_cell_count));

  dim3 block(8, 8, 4);
  for (int substep = 0; substep < substeps; ++substep) {
    for (std::size_t n = 0; n < states.size(); ++n) {
      rho_halo[n].copy_all_from(states[n].rho_d);
    }
    comm::HaloExchange::exchange_scalar(rho_halo, layout);
    compute_dry_pressure_fields(rho_halo, rho_theta_halo, pressure_fields);

    GWM_CUDA_CHECK(cudaMemset(rho_ref, 0, sizeof(real) * metrics.nz));
    GWM_CUDA_CHECK(cudaMemset(pressure_ref, 0, sizeof(real) * metrics.nz));

    for (std::size_t n = 0; n < states.size(); ++n) {
      dim3 cell_grid((states[n].rho_d.nx() + block.x - 1) / block.x,
                     (states[n].rho_d.ny() + block.y - 1) / block.y,
                     (states[n].rho_d.nz() + block.z - 1) / block.z);
      accumulate_level_reference_kernel<<<cell_grid, block>>>(
          rho_halo[n].data(), pressure_fields[n].data(), rho_ref, pressure_ref,
          states[n].rho_d.nx(), states[n].rho_d.ny(), states[n].rho_d.nz(),
          states[n].rho_d.halo());
      GWM_CUDA_CHECK(cudaGetLastError());
    }

    constexpr int level_block = 256;
    const int level_grid = (metrics.nz + level_block - 1) / level_block;
    finalize_level_reference_kernel<<<level_grid, level_block>>>(
        rho_ref, pressure_ref, metrics.nz, inv_horizontal_count);
    GWM_CUDA_CHECK(cudaGetLastError());
    sync_after_launches();

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

      dim3 xgrid((states[n].rho_d.nx() + 1 + block.x - 1) / block.x,
                 (states[n].rho_d.ny() + block.y - 1) / block.y,
                 (states[n].rho_d.nz() + block.z - 1) / block.z);
      update_fast_x_momentum_kernel<<<xgrid, block>>>(
          pressure_fields[n].data(), states[n].mom_u.storage().data(),
          states[n].rho_d.nx(), states[n].rho_d.ny(), states[n].rho_d.nz(),
          states[n].rho_d.halo(), static_cast<real>(metrics.dx), dt_fast,
          neighbors.west < 0, neighbors.east < 0);
      GWM_CUDA_CHECK(cudaGetLastError());

      dim3 ygrid((states[n].rho_d.nx() + block.x - 1) / block.x,
                 (states[n].rho_d.ny() + 1 + block.y - 1) / block.y,
                 (states[n].rho_d.nz() + block.z - 1) / block.z);
      update_fast_y_momentum_kernel<<<ygrid, block>>>(
          pressure_fields[n].data(), states[n].mom_v.storage().data(),
          states[n].rho_d.nx(), states[n].rho_d.ny(), states[n].rho_d.nz(),
          states[n].rho_d.halo(), static_cast<real>(metrics.dy), dt_fast,
          neighbors.south < 0, neighbors.north < 0);
      GWM_CUDA_CHECK(cudaGetLastError());

      dim3 zgrid((states[n].rho_d.nx() + block.x - 1) / block.x,
                 (states[n].rho_d.ny() + block.y - 1) / block.y,
                 (states[n].rho_d.nz() + 1 + block.z - 1) / block.z);
      update_fast_vertical_momentum_kernel<<<zgrid, block>>>(
          pressure_fields[n].data(), rho_halo[n].data(), pressure_ref, rho_ref,
          states[n].mom_w.storage().data(), states[n].rho_d.nx(),
          states[n].rho_d.ny(), states[n].rho_d.nz(), states[n].rho_d.halo(),
          dz, dt_fast, config.gravity);
      GWM_CUDA_CHECK(cudaGetLastError());
    }
    sync_after_launches();

    comm::HaloExchange::exchange_face(u_faces, layout);
    comm::HaloExchange::exchange_face(v_faces, layout);
    comm::HaloExchange::exchange_face(w_faces, layout);

    for (std::size_t n = 0; n < states.size(); ++n) {
      dim3 cell_grid((states[n].rho_d.nx() + block.x - 1) / block.x,
                     (states[n].rho_d.ny() + block.y - 1) / block.y,
                     (states[n].rho_d.nz() + block.z - 1) / block.z);
      update_density_from_divergence_kernel<<<cell_grid, block>>>(
          states[n].rho_d.data(), states[n].mom_u.storage().data(),
          states[n].mom_v.storage().data(), states[n].mom_w.storage().data(),
          states[n].rho_d.nx(), states[n].rho_d.ny(), states[n].rho_d.nz(),
          states[n].rho_d.halo(), static_cast<real>(metrics.dx),
          static_cast<real>(metrics.dy), dz, dt_fast);
      GWM_CUDA_CHECK(cudaGetLastError());
    }
    sync_after_launches();
  }

  GWM_CUDA_CHECK(cudaFree(rho_ref));
  GWM_CUDA_CHECK(cudaFree(pressure_ref));
}

}  // namespace gwm::dycore
