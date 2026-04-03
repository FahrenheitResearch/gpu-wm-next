#include "gwm/dycore/dry_pressure_gradient.hpp"

#include "gwm/core/cuda_utils.hpp"
#include "gwm/core/dry_thermo.hpp"
#include "gwm/comm/cartesian_topology.hpp"

namespace gwm::dycore {

namespace {

__global__ void compute_pressure_kernel(const real* rho_d, const real* rho_theta_m,
                                        real* pressure, int total_size) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= total_size) {
    return;
  }
  pressure[idx] =
      gwm::core::dry_pressure_from_rho_theta_m(rho_d[idx], rho_theta_m[idx]);
}

__global__ void compute_x_pressure_gradient_kernel(
    const real* pressure, real* mom_u_tendency, int nx, int ny, int nz, int halo,
    real dx, bool open_west, bool open_east) {
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
    mom_u_tendency[face_idx(i_face, j, k)] = 0.0f;
    return;
  }

  const real p_left = pressure[cell_idx(i_face - 1, j, k)];
  const real p_right = pressure[cell_idx(i_face, j, k)];
  mom_u_tendency[face_idx(i_face, j, k)] = -(p_right - p_left) / dx;
}

__global__ void compute_y_pressure_gradient_kernel(
    const real* pressure, real* mom_v_tendency, int nx, int ny, int nz, int halo,
    real dy, bool open_south, bool open_north) {
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
    mom_v_tendency[face_idx(i, j_face, k)] = 0.0f;
    return;
  }

  const real p_south = pressure[cell_idx(i, j_face - 1, k)];
  const real p_north = pressure[cell_idx(i, j_face, k)];
  mom_v_tendency[face_idx(i, j_face, k)] = -(p_north - p_south) / dy;
}

}  // namespace

void compute_dry_pressure_fields(
    const std::vector<state::Field3D<real>>& rho_fields,
    const std::vector<state::Field3D<real>>& rho_theta_fields,
    std::vector<state::Field3D<real>>& pressure_fields) {
  gwm::require(rho_fields.size() == rho_theta_fields.size(),
               "rho/rho_theta size mismatch in compute_dry_pressure_fields");

  if (pressure_fields.size() != rho_fields.size()) {
    pressure_fields.clear();
    pressure_fields.reserve(rho_fields.size());
    for (const auto& rho : rho_fields) {
      pressure_fields.push_back(rho.clone_empty_like("_pressure"));
    }
  }

  constexpr int block_size = 256;
  for (std::size_t n = 0; n < rho_fields.size(); ++n) {
    const int total = static_cast<int>(rho_fields[n].total_size());
    const int grid = (total + block_size - 1) / block_size;
    compute_pressure_kernel<<<grid, block_size>>>(
        rho_fields[n].data(), rho_theta_fields[n].data(),
        pressure_fields[n].data(), total);
    GWM_CUDA_CHECK(cudaGetLastError());
  }
  GWM_CUDA_CHECK(cudaDeviceSynchronize());
}

void add_horizontal_pressure_gradient_tendencies(
    const std::vector<state::Field3D<real>>& pressure_fields,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, std::vector<DrySlowTendencies>& out) {
  gwm::require(pressure_fields.size() == layout.size() &&
                   out.size() == layout.size(),
               "size mismatch in add_horizontal_pressure_gradient_tendencies");

  dim3 block(8, 8, 4);
  for (std::size_t n = 0; n < layout.size(); ++n) {
    const auto neighbors = comm::find_cartesian_neighbors(layout, layout[n]);
    dim3 xgrid((layout[n].nx_local() + 1 + block.x - 1) / block.x,
               (layout[n].ny_local() + block.y - 1) / block.y,
               (layout[n].nz + block.z - 1) / block.z);
    compute_x_pressure_gradient_kernel<<<xgrid, block>>>(
        pressure_fields[n].data(), out[n].mom_u.storage().data(),
        layout[n].nx_local(), layout[n].ny_local(), layout[n].nz, layout[n].halo,
        static_cast<real>(metrics.dx), neighbors.west < 0, neighbors.east < 0);
    GWM_CUDA_CHECK(cudaGetLastError());

    dim3 ygrid((layout[n].nx_local() + block.x - 1) / block.x,
               (layout[n].ny_local() + 1 + block.y - 1) / block.y,
               (layout[n].nz + block.z - 1) / block.z);
    compute_y_pressure_gradient_kernel<<<ygrid, block>>>(
        pressure_fields[n].data(), out[n].mom_v.storage().data(),
        layout[n].nx_local(), layout[n].ny_local(), layout[n].nz, layout[n].halo,
        static_cast<real>(metrics.dy), neighbors.south < 0, neighbors.north < 0);
    GWM_CUDA_CHECK(cudaGetLastError());
  }
  GWM_CUDA_CHECK(cudaDeviceSynchronize());
}

}  // namespace gwm::dycore
