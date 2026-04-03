#include "gwm/dycore/dry_momentum_flux.hpp"

#include "gwm/comm/cartesian_topology.hpp"
#include "gwm/core/cuda_utils.hpp"

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

__global__ void compute_cell_velocity_kernel(
    const real* rho_d, const real* mom_u, const real* mom_v, const real* mom_w,
    real* u_cell, real* v_cell, real* w_cell, int nx, int ny, int nz, int halo) {
  const int i_raw = blockIdx.x * blockDim.x + threadIdx.x;
  const int j_raw = blockIdx.y * blockDim.y + threadIdx.y;
  const int k = blockIdx.z * blockDim.z + threadIdx.z;
  const int nx_ext = nx + 2 * halo;
  const int ny_ext = ny + 2 * halo;
  if (i_raw >= nx_ext || j_raw >= ny_ext || k >= nz) {
    return;
  }

  const int i = i_raw - halo;
  const int j = j_raw - halo;

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

  const real rho = fmaxf(rho_d[cell_idx(i, j, k)], 1.0e-6f);
  u_cell[cell_idx(i, j, k)] =
      0.5f * (mom_u[xface_idx(i, j, k)] + mom_u[xface_idx(i + 1, j, k)]) / rho;
  v_cell[cell_idx(i, j, k)] =
      0.5f * (mom_v[yface_idx(i, j, k)] + mom_v[yface_idx(i, j + 1, k)]) / rho;
  w_cell[cell_idx(i, j, k)] =
      0.5f * (mom_w[zface_idx(i, j, k)] + mom_w[zface_idx(i, j, k + 1)]) / rho;
}

__device__ real upwind_flux(real adv_speed, real q_left, real q_right) {
  return adv_speed >= 0.0f ? adv_speed * q_left : adv_speed * q_right;
}

__global__ void compute_u_momentum_flux_kernel(
    const real* u_face, const real* v_face, const real* w_face, const real* u_cell,
    const real* v_cell, const real* w_cell, real* mom_u_tendency, int nx, int ny,
    int nz, int halo, const real* inv_dz_cell_metric, int metrics_nx,
    int metrics_ny, bool periodic_x, bool periodic_y, int i_begin, int j_begin,
    real dx, real dy, bool open_west, bool open_east, bool open_south,
    bool open_north) {
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
  const auto xface_idx = [&](int ii, int jj, int kk) {
    return face_storage_index(ii, jj, kk, nx + 1, ny, halo);
  };

  if ((i_face == 0 && open_west) || (i_face == nx && open_east)) {
    mom_u_tendency[xface_idx(i_face, j, k)] = 0.0f;
    return;
  }

  const real q_w = u_face[xface_idx(i_face - 1, j, k)];
  const real q_c = u_face[xface_idx(i_face, j, k)];
  const real q_e = u_face[xface_idx(i_face + 1, j, k)];
  const real q_s = u_face[xface_idx(i_face, j - 1, k)];
  const real q_n = u_face[xface_idx(i_face, j + 1, k)];
  const real q_b = (k > 0) ? u_face[xface_idx(i_face, j, k - 1)] : q_c;
  const real q_t = (k + 1 < nz) ? u_face[xface_idx(i_face, j, k + 1)] : q_c;

  const real flux_w = upwind_flux(u_cell[cell_idx(i_face - 1, j, k)], q_w, q_c);
  const real flux_e = upwind_flux(u_cell[cell_idx(i_face, j, k)], q_c, q_e);

  real flux_s = 0.0f;
  real flux_n = 0.0f;
  if (!(j == 0 && open_south)) {
    const real v_s = 0.5f * (v_cell[cell_idx(i_face - 1, j - 1, k)] +
                             v_cell[cell_idx(i_face, j - 1, k)]);
    flux_s = upwind_flux(v_s, q_s, q_c);
  }
  if (!(j + 1 == ny && open_north)) {
    const real v_n = 0.5f * (v_cell[cell_idx(i_face - 1, j, k)] +
                             v_cell[cell_idx(i_face, j, k)]);
    flux_n = upwind_flux(v_n, q_c, q_n);
  }

  real flux_b = 0.0f;
  real flux_t = 0.0f;
  if (k > 0) {
    const real w_b = 0.5f * (w_cell[cell_idx(i_face - 1, j, k - 1)] +
                             w_cell[cell_idx(i_face, j, k - 1)]);
    flux_b = upwind_flux(w_b, q_b, q_c);
  }
  if (k + 1 < nz) {
    const real w_t = 0.5f * (w_cell[cell_idx(i_face - 1, j, k)] +
                             w_cell[cell_idx(i_face, j, k)]);
    flux_t = upwind_flux(w_t, q_c, q_t);
  }

  const int gj = wrap_metric_index(j_begin + j, metrics_ny, periodic_y);
  const int gi_left =
      wrap_metric_index(i_begin + i_face - 1, metrics_nx, periodic_x);
  const int gi_right =
      wrap_metric_index(i_begin + i_face, metrics_nx, periodic_x);
  const real inv_dz =
      0.5f *
      (inv_dz_cell_metric[metric_storage_index(gi_left, gj, k, metrics_nx, metrics_ny)] +
       inv_dz_cell_metric[metric_storage_index(gi_right, gj, k, metrics_nx, metrics_ny)]);

  mom_u_tendency[xface_idx(i_face, j, k)] =
      -((flux_e - flux_w) / dx + (flux_n - flux_s) / dy +
        (flux_t - flux_b) * inv_dz);
}

__global__ void compute_v_momentum_flux_kernel(
    const real* u_face, const real* v_face, const real* w_face, const real* u_cell,
    const real* v_cell, const real* w_cell, real* mom_v_tendency, int nx, int ny,
    int nz, int halo, const real* inv_dz_cell_metric, int metrics_nx,
    int metrics_ny, bool periodic_x, bool periodic_y, int i_begin, int j_begin,
    real dx, real dy, bool open_west, bool open_east, bool open_south,
    bool open_north) {
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
  const auto yface_idx = [&](int ii, int jj, int kk) {
    return face_storage_index(ii, jj, kk, nx, ny + 1, halo);
  };

  if ((j_face == 0 && open_south) || (j_face == ny && open_north)) {
    mom_v_tendency[yface_idx(i, j_face, k)] = 0.0f;
    return;
  }

  const real q_w = v_face[yface_idx(i - 1, j_face, k)];
  const real q_c = v_face[yface_idx(i, j_face, k)];
  const real q_e = v_face[yface_idx(i + 1, j_face, k)];
  const real q_s = v_face[yface_idx(i, j_face - 1, k)];
  const real q_n = v_face[yface_idx(i, j_face + 1, k)];
  const real q_b = (k > 0) ? v_face[yface_idx(i, j_face, k - 1)] : q_c;
  const real q_t = (k + 1 < nz) ? v_face[yface_idx(i, j_face, k + 1)] : q_c;

  real flux_w = 0.0f;
  real flux_e = 0.0f;
  if (!(i == 0 && open_west)) {
    const real u_w = 0.5f * (u_cell[cell_idx(i - 1, j_face - 1, k)] +
                             u_cell[cell_idx(i - 1, j_face, k)]);
    flux_w = upwind_flux(u_w, q_w, q_c);
  }
  if (!(i + 1 == nx && open_east)) {
    const real u_e = 0.5f * (u_cell[cell_idx(i, j_face - 1, k)] +
                             u_cell[cell_idx(i, j_face, k)]);
    flux_e = upwind_flux(u_e, q_c, q_e);
  }

  const real flux_s = upwind_flux(v_cell[cell_idx(i, j_face - 1, k)], q_s, q_c);
  const real flux_n = upwind_flux(v_cell[cell_idx(i, j_face, k)], q_c, q_n);

  real flux_b = 0.0f;
  real flux_t = 0.0f;
  if (k > 0) {
    const real w_b = 0.5f * (w_cell[cell_idx(i, j_face - 1, k - 1)] +
                             w_cell[cell_idx(i, j_face, k - 1)]);
    flux_b = upwind_flux(w_b, q_b, q_c);
  }
  if (k + 1 < nz) {
    const real w_t = 0.5f * (w_cell[cell_idx(i, j_face - 1, k)] +
                             w_cell[cell_idx(i, j_face, k)]);
    flux_t = upwind_flux(w_t, q_c, q_t);
  }

  const int gi = wrap_metric_index(i_begin + i, metrics_nx, periodic_x);
  const int gj_south =
      wrap_metric_index(j_begin + j_face - 1, metrics_ny, periodic_y);
  const int gj_north =
      wrap_metric_index(j_begin + j_face, metrics_ny, periodic_y);
  const real inv_dz =
      0.5f *
      (inv_dz_cell_metric[metric_storage_index(gi, gj_south, k, metrics_nx, metrics_ny)] +
       inv_dz_cell_metric[metric_storage_index(gi, gj_north, k, metrics_nx, metrics_ny)]);

  mom_v_tendency[yface_idx(i, j_face, k)] =
      -((flux_e - flux_w) / dx + (flux_n - flux_s) / dy +
        (flux_t - flux_b) * inv_dz);
}

__global__ void compute_w_momentum_flux_kernel(
    const real* u_face, const real* v_face, const real* w_face, const real* u_cell,
    const real* v_cell, const real* w_cell, real* mom_w_tendency, int nx, int ny,
    int nz, int halo, const real* inv_dz_face_metric, int metrics_nx,
    int metrics_ny, bool periodic_x, bool periodic_y, int i_begin, int j_begin,
    real dx, real dy, bool open_west, bool open_east, bool open_south,
    bool open_north) {
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
    mom_w_tendency[zface_idx(i, j, k_face)] = 0.0f;
    return;
  }

  const int k_lower = k_face - 1;
  const int k_upper = k_face;

  const real q_w = w_face[zface_idx(i - 1, j, k_face)];
  const real q_c = w_face[zface_idx(i, j, k_face)];
  const real q_e = w_face[zface_idx(i + 1, j, k_face)];
  const real q_s = w_face[zface_idx(i, j - 1, k_face)];
  const real q_n = w_face[zface_idx(i, j + 1, k_face)];
  const real q_b = w_face[zface_idx(i, j, k_face - 1)];
  const real q_t = w_face[zface_idx(i, j, k_face + 1)];

  real flux_w = 0.0f;
  real flux_e = 0.0f;
  if (!(i == 0 && open_west)) {
    const real u_w = 0.5f * (u_cell[cell_idx(i - 1, j, k_lower)] +
                             u_cell[cell_idx(i - 1, j, k_upper)]);
    flux_w = upwind_flux(u_w, q_w, q_c);
  }
  if (!(i + 1 == nx && open_east)) {
    const real u_e = 0.5f * (u_cell[cell_idx(i, j, k_lower)] +
                             u_cell[cell_idx(i, j, k_upper)]);
    flux_e = upwind_flux(u_e, q_c, q_e);
  }

  real flux_s = 0.0f;
  real flux_n = 0.0f;
  if (!(j == 0 && open_south)) {
    const real v_s = 0.5f * (v_cell[cell_idx(i, j - 1, k_lower)] +
                             v_cell[cell_idx(i, j - 1, k_upper)]);
    flux_s = upwind_flux(v_s, q_s, q_c);
  }
  if (!(j + 1 == ny && open_north)) {
    const real v_n = 0.5f * (v_cell[cell_idx(i, j, k_lower)] +
                             v_cell[cell_idx(i, j, k_upper)]);
    flux_n = upwind_flux(v_n, q_c, q_n);
  }

  const real flux_b = upwind_flux(w_cell[cell_idx(i, j, k_lower)], q_b, q_c);
  const real flux_t = upwind_flux(w_cell[cell_idx(i, j, k_upper)], q_c, q_t);
  const real inv_dz =
      inv_dz_face_metric[metric_storage_index(
          wrap_metric_index(i_begin + i, metrics_nx, periodic_x),
          wrap_metric_index(j_begin + j, metrics_ny, periodic_y), k_face,
          metrics_nx, metrics_ny)];

  mom_w_tendency[zface_idx(i, j, k_face)] =
      -((flux_e - flux_w) / dx + (flux_n - flux_s) / dy +
        (flux_t - flux_b) * inv_dz);
}

void sync_after_launches() { GWM_CUDA_CHECK(cudaDeviceSynchronize()); }

}  // namespace

void add_dry_momentum_flux_tendencies(
    const std::vector<state::Field3D<real>>& rho_fields,
    const std::vector<state::FaceField<real>>& mom_u_fields,
    const std::vector<state::FaceField<real>>& mom_v_fields,
    const std::vector<state::FaceField<real>>& mom_w_fields,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, std::vector<DrySlowTendencies>& out) {
  gwm::require(rho_fields.size() == layout.size() &&
                   mom_u_fields.size() == layout.size() &&
                   mom_v_fields.size() == layout.size() &&
                   mom_w_fields.size() == layout.size() &&
                   out.size() == layout.size(),
               "size mismatch in add_dry_momentum_flux_tendencies");

  std::vector<state::Field3D<real>> u_cell;
  std::vector<state::Field3D<real>> v_cell;
  std::vector<state::Field3D<real>> w_cell;
  u_cell.reserve(layout.size());
  v_cell.reserve(layout.size());
  w_cell.reserve(layout.size());
  for (const auto& rho : rho_fields) {
    u_cell.push_back(rho.clone_empty_like("_u_cell"));
    v_cell.push_back(rho.clone_empty_like("_v_cell"));
    w_cell.push_back(rho.clone_empty_like("_w_cell"));
  }

  dim3 block(8, 8, 4);
  for (std::size_t n = 0; n < layout.size(); ++n) {
    const int nx_ext = rho_fields[n].nx() + 2 * rho_fields[n].halo();
    const int ny_ext = rho_fields[n].ny() + 2 * rho_fields[n].halo();
    dim3 vel_grid((nx_ext + block.x - 1) / block.x,
                  (ny_ext + block.y - 1) / block.y,
                  (rho_fields[n].nz() + block.z - 1) / block.z);
    compute_cell_velocity_kernel<<<vel_grid, block>>>(
        rho_fields[n].data(), mom_u_fields[n].storage().data(),
        mom_v_fields[n].storage().data(), mom_w_fields[n].storage().data(),
        u_cell[n].data(), v_cell[n].data(), w_cell[n].data(), rho_fields[n].nx(),
        rho_fields[n].ny(), rho_fields[n].nz(), rho_fields[n].halo());
    GWM_CUDA_CHECK(cudaGetLastError());
  }
  sync_after_launches();

  for (std::size_t n = 0; n < layout.size(); ++n) {
    const auto neighbors = comm::find_cartesian_neighbors(layout, layout[n]);

    dim3 u_grid((layout[n].nx_local() + 1 + block.x - 1) / block.x,
                (layout[n].ny_local() + block.y - 1) / block.y,
                (layout[n].nz + block.z - 1) / block.z);
    compute_u_momentum_flux_kernel<<<u_grid, block>>>(
        mom_u_fields[n].storage().data(), mom_v_fields[n].storage().data(),
        mom_w_fields[n].storage().data(), u_cell[n].data(), v_cell[n].data(),
        w_cell[n].data(), out[n].mom_u.storage().data(), layout[n].nx_local(),
        layout[n].ny_local(), layout[n].nz, layout[n].halo,
        metrics.inv_dz_cell_metric.data(), metrics.nx, metrics.ny,
        metrics.periodic_x, metrics.periodic_y, layout[n].i_begin,
        layout[n].j_begin, static_cast<real>(metrics.dx),
        static_cast<real>(metrics.dy), neighbors.west < 0, neighbors.east < 0,
        neighbors.south < 0, neighbors.north < 0);
    GWM_CUDA_CHECK(cudaGetLastError());

    dim3 v_grid((layout[n].nx_local() + block.x - 1) / block.x,
                (layout[n].ny_local() + 1 + block.y - 1) / block.y,
                (layout[n].nz + block.z - 1) / block.z);
    compute_v_momentum_flux_kernel<<<v_grid, block>>>(
        mom_u_fields[n].storage().data(), mom_v_fields[n].storage().data(),
        mom_w_fields[n].storage().data(), u_cell[n].data(), v_cell[n].data(),
        w_cell[n].data(), out[n].mom_v.storage().data(), layout[n].nx_local(),
        layout[n].ny_local(), layout[n].nz, layout[n].halo,
        metrics.inv_dz_cell_metric.data(), metrics.nx, metrics.ny,
        metrics.periodic_x, metrics.periodic_y, layout[n].i_begin,
        layout[n].j_begin, static_cast<real>(metrics.dx),
        static_cast<real>(metrics.dy), neighbors.west < 0, neighbors.east < 0,
        neighbors.south < 0, neighbors.north < 0);
    GWM_CUDA_CHECK(cudaGetLastError());

    dim3 w_grid((layout[n].nx_local() + block.x - 1) / block.x,
                (layout[n].ny_local() + block.y - 1) / block.y,
                (layout[n].nz + 1 + block.z - 1) / block.z);
    compute_w_momentum_flux_kernel<<<w_grid, block>>>(
        mom_u_fields[n].storage().data(), mom_v_fields[n].storage().data(),
        mom_w_fields[n].storage().data(), u_cell[n].data(), v_cell[n].data(),
        w_cell[n].data(), out[n].mom_w.storage().data(), layout[n].nx_local(),
        layout[n].ny_local(), layout[n].nz, layout[n].halo,
        metrics.inv_dz_face_metric.data(), metrics.nx, metrics.ny,
        metrics.periodic_x, metrics.periodic_y, layout[n].i_begin,
        layout[n].j_begin, static_cast<real>(metrics.dx),
        static_cast<real>(metrics.dy), neighbors.west < 0, neighbors.east < 0,
        neighbors.south < 0, neighbors.north < 0);
    GWM_CUDA_CHECK(cudaGetLastError());
  }
  sync_after_launches();
}

}  // namespace gwm::dycore
