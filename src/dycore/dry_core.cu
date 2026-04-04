#include "gwm/dycore/dry_core.hpp"

#include <algorithm>
#include <cmath>

#include "gwm/comm/collectives.hpp"
#include "gwm/comm/halo_exchange.hpp"
#include "gwm/core/cuda_utils.hpp"
#include "gwm/core/dry_thermo.hpp"
#include "gwm/dycore/dry_fast_modes.hpp"
#include "gwm/dycore/dry_momentum_flux.hpp"
#include "gwm/dycore/dry_pressure_gradient.hpp"

namespace gwm::dycore {

namespace {

constexpr real kGravity = 9.81f;

real hydrostatic_rho_at_height(real rho_surface, real theta_ref, real z,
                               real gravity) {
  const real p_surface = gwm::core::dry_pressure_from_rho_theta_m(
      rho_surface, rho_surface * theta_ref);
  const real exner_surface =
      powf(fmaxf(p_surface / gwm::core::kReferencePressure, 1.0e-6f),
           gwm::core::kKappa);
  const real exner =
      fmaxf(1.0e-4f, exner_surface -
                         gravity * z / (gwm::core::dry_cp() * theta_ref));
  const real pressure =
      gwm::core::kReferencePressure * powf(exner, 1.0f / gwm::core::kKappa);
  return gwm::core::dry_rho_from_pressure_theta(pressure, theta_ref);
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

__host__ __device__ std::size_t metric_storage_index(int i, int j, int k,
                                                     int nx, int ny) {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

__global__ void fill_tendencies_zero_kernel(real* data, int total_size) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < total_size) {
    data[idx] = 0.0f;
  }
}

__global__ void compute_cell_tendencies_kernel(
    const real* rho_d, const real* rho_theta_m, const real* mom_u,
    const real* mom_v, const real* mom_w, real* rho_tendency,
    real* rho_theta_tendency, const real* inv_dz_cell_metric, int metrics_nx,
    int metrics_ny, int i_begin, int j_begin, int nx, int ny, int nz,
    int halo, real dx, real dy) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  const int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (i >= nx || j >= ny || k >= nz) {
    return;
  }

  const int cell_pitch_x = nx + 2 * halo;
  const int cell_pitch_y = ny + 2 * halo;
  const auto cell_idx = [&](int ii, int jj, int kk) {
    return (kk * cell_pitch_y + (jj + halo)) * cell_pitch_x + (ii + halo);
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

  const real mass_flux_w = mom_u[xface_idx(i, j, k)];
  const real mass_flux_e = mom_u[xface_idx(i + 1, j, k)];
  const real mass_flux_s = mom_v[yface_idx(i, j, k)];
  const real mass_flux_n = mom_v[yface_idx(i, j + 1, k)];
  const real mass_flux_b = mom_w[zface_idx(i, j, k)];
  const real mass_flux_t = mom_w[zface_idx(i, j, k + 1)];
  const real inv_dz =
      inv_dz_cell_metric[metric_storage_index(i_begin + i, j_begin + j, k,
                                              metrics_nx, metrics_ny)];

  const real div_mass = (mass_flux_e - mass_flux_w) / dx +
                        (mass_flux_n - mass_flux_s) / dy +
                        (mass_flux_t - mass_flux_b) * inv_dz;
  rho_tendency[cell_idx(i, j, k)] = -div_mass;

  const real rho_c = rho_d[cell_idx(i, j, k)];
  const real theta_c =
      gwm::core::dry_theta_from_conserved(rho_c, rho_theta_m[cell_idx(i, j, k)]);

  const real theta_w =
      gwm::core::dry_theta_from_conserved(rho_d[cell_idx(i - 1, j, k)],
                                          rho_theta_m[cell_idx(i - 1, j, k)]);
  const real theta_e =
      gwm::core::dry_theta_from_conserved(rho_d[cell_idx(i + 1, j, k)],
                                          rho_theta_m[cell_idx(i + 1, j, k)]);
  const real theta_s =
      gwm::core::dry_theta_from_conserved(rho_d[cell_idx(i, j - 1, k)],
                                          rho_theta_m[cell_idx(i, j - 1, k)]);
  const real theta_n =
      gwm::core::dry_theta_from_conserved(rho_d[cell_idx(i, j + 1, k)],
                                          rho_theta_m[cell_idx(i, j + 1, k)]);

  const real theta_b =
      (k > 0) ? gwm::core::dry_theta_from_conserved(
                    rho_d[cell_idx(i, j, k - 1)],
                    rho_theta_m[cell_idx(i, j, k - 1)])
              : theta_c;
  const real theta_t =
      (k + 1 < nz) ? gwm::core::dry_theta_from_conserved(
                         rho_d[cell_idx(i, j, k + 1)],
                         rho_theta_m[cell_idx(i, j, k + 1)])
                   : theta_c;

  const real flux_theta_w = mass_flux_w * (mass_flux_w >= 0.0f ? theta_w : theta_c);
  const real flux_theta_e = mass_flux_e * (mass_flux_e >= 0.0f ? theta_c : theta_e);
  const real flux_theta_s = mass_flux_s * (mass_flux_s >= 0.0f ? theta_s : theta_c);
  const real flux_theta_n = mass_flux_n * (mass_flux_n >= 0.0f ? theta_c : theta_n);
  const real flux_theta_b = mass_flux_b * (mass_flux_b >= 0.0f ? theta_b : theta_c);
  const real flux_theta_t = mass_flux_t * (mass_flux_t >= 0.0f ? theta_c : theta_t);

  const real div_rho_theta = (flux_theta_e - flux_theta_w) / dx +
                             (flux_theta_n - flux_theta_s) / dy +
                             (flux_theta_t - flux_theta_b) * inv_dz;
  rho_theta_tendency[cell_idx(i, j, k)] = -div_rho_theta;
}

__global__ void compute_vertical_buoyancy_kernel(
    const real* rho_d, const real* rho_theta_m, const real* theta_ref,
    real* mom_w_tendency, int nx, int ny, int nz, int halo, real gravity) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  const int k_face = blockIdx.z * blockDim.z + threadIdx.z;
  if (i >= nx || j >= ny || k_face > nz) {
    return;
  }

  const int cell_pitch_x = nx + 2 * halo;
  const int cell_pitch_y = ny + 2 * halo;
  const auto cell_idx = [&](int ii, int jj, int kk) {
    return (kk * cell_pitch_y + (jj + halo)) * cell_pitch_x + (ii + halo);
  };
  const auto zface_idx = [&](int ii, int jj, int kk) {
    return face_storage_index(ii, jj, kk, nx, ny, halo);
  };

  if (k_face == 0 || k_face == nz) {
    return;
  }

  const int k_lower = k_face - 1;
  const int k_upper = k_face;

  const real rho_lower = rho_d[cell_idx(i, j, k_lower)];
  const real rho_upper = rho_d[cell_idx(i, j, k_upper)];
  const real theta_lower =
      rho_theta_m[cell_idx(i, j, k_lower)] / fmaxf(rho_lower, 1.0e-6f);
  const real theta_upper =
      rho_theta_m[cell_idx(i, j, k_upper)] / fmaxf(rho_upper, 1.0e-6f);

  const real theta_face = 0.5f * (theta_lower + theta_upper);
  const real theta_ref_face =
      0.5f * (theta_ref[k_lower] + theta_ref[k_upper]);
  const real rho_face = 0.5f * (rho_lower + rho_upper);
  const real buoyancy =
      gravity * (theta_face - theta_ref_face) / fmaxf(theta_ref_face, 1.0e-6f);

  mom_w_tendency[zface_idx(i, j, k_face)] += rho_face * buoyancy;
}

__global__ void update_cell_field_kernel(real* out, const real* base,
                                         const real* tendency, int total_size,
                                         real scale, real clip_floor) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= total_size) {
    return;
  }
  out[idx] = fmaxf(clip_floor, base[idx] + scale * tendency[idx]);
}

__global__ void blend_cell_field_kernel(real* out, const real* a, const real* b,
                                        int total_size, real wa, real wb) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < total_size) {
    out[idx] = wa * a[idx] + wb * b[idx];
  }
}

void launch_zero_fill(state::Field3D<real>& field) {
  constexpr int block_size = 256;
  const int total = static_cast<int>(field.total_size());
  const int grid = (total + block_size - 1) / block_size;
  fill_tendencies_zero_kernel<<<grid, block_size>>>(field.data(), total);
  GWM_CUDA_CHECK(cudaGetLastError());
}

void launch_zero_fill(state::FaceField<real>& field) {
  launch_zero_fill(field.storage());
}

void update_stage_field(state::Field3D<real>& out, const state::Field3D<real>& base,
                        const state::Field3D<real>& tendency, real dt,
                        real clip_floor) {
  constexpr int block_size = 256;
  const int total = static_cast<int>(out.total_size());
  const int grid = (total + block_size - 1) / block_size;
  update_cell_field_kernel<<<grid, block_size>>>(out.data(), base.data(),
                                                 tendency.data(), total, dt,
                                                 clip_floor);
  GWM_CUDA_CHECK(cudaGetLastError());
}

void update_stage_field(state::FaceField<real>& out,
                        const state::FaceField<real>& base,
                        const state::FaceField<real>& tendency, real dt,
                        real clip_floor) {
  update_stage_field(out.storage(), base.storage(), tendency.storage(), dt,
                     clip_floor);
}

void blend_stage_field(state::Field3D<real>& out, const state::Field3D<real>& a,
                       const state::Field3D<real>& b, real wa, real wb) {
  constexpr int block_size = 256;
  const int total = static_cast<int>(out.total_size());
  const int grid = (total + block_size - 1) / block_size;
  blend_cell_field_kernel<<<grid, block_size>>>(out.data(), a.data(), b.data(),
                                                total, wa, wb);
  GWM_CUDA_CHECK(cudaGetLastError());
}

void blend_stage_field(state::FaceField<real>& out,
                       const state::FaceField<real>& a,
                       const state::FaceField<real>& b, real wa, real wb) {
  blend_stage_field(out.storage(), a.storage(), b.storage(), wa, wb);
}

void sync_after_launches() { GWM_CUDA_CHECK(cudaDeviceSynchronize()); }

}  // namespace

DryState::DryState(int nx, int ny, int nz, int halo,
                   const std::string& label_prefix)
    : rho_d(nx, ny, nz, halo, label_prefix + "_rho_d"),
      rho_theta_m(nx, ny, nz, halo, label_prefix + "_rho_theta_m"),
      mom_u(nx, ny, nz, halo, state::FaceOrientation::X,
            label_prefix + "_mom_u"),
      mom_v(nx, ny, nz, halo, state::FaceOrientation::Y,
            label_prefix + "_mom_v"),
      mom_w(nx, ny, nz, halo, state::FaceOrientation::Z,
            label_prefix + "_mom_w") {}

DryState DryState::clone_empty_like(const std::string& suffix) const {
  return DryState(rho_d.nx(), rho_d.ny(), rho_d.nz(), rho_d.halo(),
                  suffix.empty() ? "dry_state" : suffix);
}

void DryState::copy_all_from(const DryState& other) {
  rho_d.copy_all_from(other.rho_d);
  rho_theta_m.copy_all_from(other.rho_theta_m);
  mom_u.storage().copy_all_from(other.mom_u.storage());
  mom_v.storage().copy_all_from(other.mom_v.storage());
  mom_w.storage().copy_all_from(other.mom_w.storage());
}

void DryState::fill_constant(real rho_value, real theta_value, real u_value,
                             real v_value, real w_value) {
  rho_d.fill(rho_value);
  rho_theta_m.fill(rho_value * theta_value);
  mom_u.storage().fill(rho_value * u_value);
  mom_v.storage().fill(rho_value * v_value);
  mom_w.storage().fill(rho_value * w_value);
}

double DryState::total_dry_mass() const { return rho_d.owned_sum(); }

double DryState::total_rho_theta_m() const { return rho_theta_m.owned_sum(); }

double DryState::total_vertical_momentum() const {
  return mom_w.storage().owned_sum();
}

DrySlowTendencies::DrySlowTendencies(int nx, int ny, int nz, int halo,
                                     const std::string& label_prefix)
    : rho_d(nx, ny, nz, halo, label_prefix + "_rho_d"),
      rho_theta_m(nx, ny, nz, halo, label_prefix + "_rho_theta_m"),
      mom_u(nx, ny, nz, halo, state::FaceOrientation::X,
            label_prefix + "_mom_u"),
      mom_v(nx, ny, nz, halo, state::FaceOrientation::Y,
            label_prefix + "_mom_v"),
      mom_w(nx, ny, nz, halo, state::FaceOrientation::Z,
            label_prefix + "_mom_w") {}

void DrySlowTendencies::fill_zero() {
  launch_zero_fill(rho_d);
  launch_zero_fill(rho_theta_m);
  launch_zero_fill(mom_u);
  launch_zero_fill(mom_v);
  launch_zero_fill(mom_w);
  sync_after_launches();
}

void NullBoundaryUpdater::apply(
    std::vector<DryState>& /*states*/,
    const std::vector<domain::SubdomainDescriptor>& /*layout*/,
    real /*sim_time*/) {}

void LocalSplitExplicitFastMode::apply(
    std::vector<DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const DryStepperConfig& config) {
  apply_local_split_explicit_fast_modes(states, layout, metrics, config);
}

std::vector<DryState> make_constant_dry_state(
    const std::vector<domain::SubdomainDescriptor>& layout, real rho_value,
    real theta_value, real u_value, real v_value, real w_value,
    const std::string& label_prefix) {
  std::vector<DryState> states;
  states.reserve(layout.size());
  for (const auto& desc : layout) {
    DryState state(desc.nx_local(), desc.ny_local(), desc.nz, desc.halo,
                   label_prefix + "_rank_" + std::to_string(desc.rank));
    state.fill_constant(rho_value, theta_value, u_value, v_value, w_value);
    states.push_back(std::move(state));
  }
  return states;
}

std::vector<DryState> make_hydrostatic_rest_state(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, real rho_surface, real theta_ref,
    real density_scale_height, const std::string& label_prefix) {
  (void)density_scale_height;
  gwm::require(static_cast<int>(metrics.z_centers.size()) == metrics.nz,
               "Grid metrics must contain z-centers for hydrostatic state");

  std::vector<DryState> states;
  states.reserve(layout.size());
  for (const auto& desc : layout) {
    DryState state(desc.nx_local(), desc.ny_local(), desc.nz, desc.halo,
                   label_prefix + "_rank_" + std::to_string(desc.rank));
    state.fill_constant(0.0f, theta_ref, 0.0f, 0.0f, 0.0f);

    for (int k = 0; k < desc.nz; ++k) {
      for (int j = 0; j < desc.ny_local(); ++j) {
        const int j_global = desc.j_begin + j;
        for (int i = 0; i < desc.nx_local(); ++i) {
          const int i_global = desc.i_begin + i;
          const real rho_k = hydrostatic_rho_at_height(
              rho_surface, theta_ref, metrics.z_center(i_global, j_global, k),
              kGravity);
          state.rho_d(i, j, k) = rho_k;
          state.rho_theta_m(i, j, k) = rho_k * theta_ref;
        }
      }
    }

    state.mom_u.storage().fill(0.0f);
    state.mom_v.storage().fill(0.0f);
    state.mom_w.storage().fill(0.0f);
    states.push_back(std::move(state));
  }
  return states;
}

void compute_slow_tendencies(
    const std::vector<DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const DryStepperConfig& config,
    std::vector<DrySlowTendencies>& out) {
  gwm::require(states.size() == layout.size(),
               "State/layout size mismatch in compute_slow_tendencies");

  if (out.size() != states.size()) {
    out.clear();
    out.reserve(states.size());
    for (const auto& state : states) {
      out.emplace_back(state.rho_d.nx(), state.rho_d.ny(), state.rho_d.nz(),
                       state.rho_d.halo(), "slow_tendency");
    }
  }

  std::vector<state::Field3D<real>> rho_fields;
  std::vector<state::Field3D<real>> rho_theta_fields;
  std::vector<state::Field3D<real>> pressure_fields;
  std::vector<state::FaceField<real>> mom_u_fields;
  std::vector<state::FaceField<real>> mom_v_fields;
  std::vector<state::FaceField<real>> mom_w_fields;
  rho_fields.reserve(states.size());
  rho_theta_fields.reserve(states.size());
  pressure_fields.reserve(states.size());
  mom_u_fields.reserve(states.size());
  mom_v_fields.reserve(states.size());
  mom_w_fields.reserve(states.size());
  for (const auto& state : states) {
    rho_fields.push_back(state.rho_d.clone_empty_like("_rho_halo"));
    rho_fields.back().copy_all_from(state.rho_d);
    rho_theta_fields.push_back(state.rho_theta_m.clone_empty_like("_theta_halo"));
    rho_theta_fields.back().copy_all_from(state.rho_theta_m);
    mom_u_fields.emplace_back(state.rho_d.nx(), state.rho_d.ny(), state.rho_d.nz(),
                              state.rho_d.halo(), state::FaceOrientation::X,
                              "mom_u_halo");
    mom_u_fields.back().storage().copy_all_from(state.mom_u.storage());
    mom_v_fields.emplace_back(state.rho_d.nx(), state.rho_d.ny(), state.rho_d.nz(),
                              state.rho_d.halo(), state::FaceOrientation::Y,
                              "mom_v_halo");
    mom_v_fields.back().storage().copy_all_from(state.mom_v.storage());
    mom_w_fields.emplace_back(state.rho_d.nx(), state.rho_d.ny(), state.rho_d.nz(),
                              state.rho_d.halo(), state::FaceOrientation::Z,
                              "mom_w_halo");
    mom_w_fields.back().storage().copy_all_from(state.mom_w.storage());
  }

  comm::HaloExchange::exchange_scalar(rho_fields, layout);
  comm::HaloExchange::exchange_scalar(rho_theta_fields, layout);
  comm::HaloExchange::exchange_face(mom_u_fields, layout);
  comm::HaloExchange::exchange_face(mom_v_fields, layout);
  comm::HaloExchange::exchange_face(mom_w_fields, layout);
  comm::HaloExchange::synchronize_owned_face_interfaces(mom_u_fields, layout);
  comm::HaloExchange::synchronize_owned_face_interfaces(mom_v_fields, layout);
  comm::HaloExchange::exchange_face(mom_u_fields, layout);
  comm::HaloExchange::exchange_face(mom_v_fields, layout);
  compute_dry_pressure_fields(rho_fields, rho_theta_fields, pressure_fields);

  std::vector<double> theta_sum(static_cast<std::size_t>(metrics.nz), 0.0);
  std::vector<std::uint64_t> theta_count(static_cast<std::size_t>(metrics.nz),
                                         0);
  for (const auto& state : states) {
    for (int k = 0; k < state.rho_d.nz(); ++k) {
      for (int j = 0; j < state.rho_d.ny(); ++j) {
        for (int i = 0; i < state.rho_d.nx(); ++i) {
          const real rho = state.rho_d(i, j, k);
          const real theta = state.rho_theta_m(i, j, k) / std::max(rho, 1.0e-6f);
          theta_sum[static_cast<std::size_t>(k)] += static_cast<double>(theta);
          theta_count[static_cast<std::size_t>(k)] += 1u;
        }
      }
    }
  }

  comm::allreduce_sum_in_place(theta_sum);
  comm::allreduce_sum_in_place(theta_count);

  real* theta_ref_device = nullptr;
  GWM_CUDA_CHECK(
      cudaMallocManaged(&theta_ref_device, sizeof(real) * metrics.nz));
  for (int k = 0; k < metrics.nz; ++k) {
    const auto idx = static_cast<std::size_t>(k);
    theta_ref_device[k] = static_cast<real>(
        theta_sum[idx] /
        static_cast<double>(std::max<std::uint64_t>(1u, theta_count[idx])));
  }

  for (std::size_t n = 0; n < states.size(); ++n) {
    out[n].fill_zero();

    dim3 block(8, 8, 4);
    dim3 grid((states[n].rho_d.nx() + block.x - 1) / block.x,
              (states[n].rho_d.ny() + block.y - 1) / block.y,
              (states[n].rho_d.nz() + block.z - 1) / block.z);

    compute_cell_tendencies_kernel<<<grid, block>>>(
        rho_fields[n].data(), rho_theta_fields[n].data(),
        mom_u_fields[n].storage().data(), mom_v_fields[n].storage().data(),
        mom_w_fields[n].storage().data(), out[n].rho_d.data(),
        out[n].rho_theta_m.data(), metrics.inv_dz_cell_metric.data(),
        metrics.nx, metrics.ny, layout[n].i_begin, layout[n].j_begin,
        states[n].rho_d.nx(), states[n].rho_d.ny(), states[n].rho_d.nz(),
        states[n].rho_d.halo(), static_cast<real>(metrics.dx),
        static_cast<real>(metrics.dy));
    GWM_CUDA_CHECK(cudaGetLastError());
  }

  add_dry_momentum_flux_tendencies(rho_fields, mom_u_fields, mom_v_fields,
                                   mom_w_fields, layout, metrics, out);

  add_horizontal_pressure_gradient_tendencies(pressure_fields, layout, metrics,
                                              out);

  for (std::size_t n = 0; n < states.size(); ++n) {
    dim3 block(8, 8, 4);
    dim3 face_grid((states[n].rho_d.nx() + block.x - 1) / block.x,
                   (states[n].rho_d.ny() + block.y - 1) / block.y,
                   (states[n].rho_d.nz() + 1 + block.z - 1) / block.z);
    compute_vertical_buoyancy_kernel<<<face_grid, block>>>(
        rho_fields[n].data(), rho_theta_fields[n].data(), theta_ref_device,
        out[n].mom_w.storage().data(), states[n].rho_d.nx(),
        states[n].rho_d.ny(), states[n].rho_d.nz(), states[n].rho_d.halo(),
        config.gravity);
    GWM_CUDA_CHECK(cudaGetLastError());
  }
  sync_after_launches();
  GWM_CUDA_CHECK(cudaFree(theta_ref_device));
}

void advance_dry_state_ssprk3(
    std::vector<DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const DryStepperConfig& config,
    BoundaryUpdater& boundary_updater, FastModeIntegrator& fast_modes) {
  gwm::require(states.size() == layout.size(),
               "State/layout size mismatch in advance_dry_state_ssprk3");

  std::vector<DryState> q0;
  std::vector<DryState> q1;
  std::vector<DryState> q2;
  std::vector<DryState> tmp;
  q0.reserve(states.size());
  q1.reserve(states.size());
  q2.reserve(states.size());
  tmp.reserve(states.size());
  for (const auto& state : states) {
    q0.push_back(state.clone_empty_like("q0"));
    q0.back().copy_all_from(state);
    q1.push_back(state.clone_empty_like("q1"));
    q2.push_back(state.clone_empty_like("q2"));
    tmp.push_back(state.clone_empty_like("tmp"));
  }

  std::vector<DrySlowTendencies> tendencies;

  boundary_updater.apply(states, layout, 0.0f);
  compute_slow_tendencies(states, layout, metrics, config, tendencies);
  for (std::size_t n = 0; n < states.size(); ++n) {
    q1[n].copy_all_from(states[n]);
    update_stage_field(q1[n].rho_d, states[n].rho_d, tendencies[n].rho_d,
                       config.dt, 1.0e-6f);
    update_stage_field(q1[n].rho_theta_m, states[n].rho_theta_m,
                       tendencies[n].rho_theta_m, config.dt, 1.0e-6f);
    update_stage_field(q1[n].mom_u, states[n].mom_u, tendencies[n].mom_u,
                       config.dt, -1.0e30f);
    update_stage_field(q1[n].mom_v, states[n].mom_v, tendencies[n].mom_v,
                       config.dt, -1.0e30f);
    update_stage_field(q1[n].mom_w, states[n].mom_w, tendencies[n].mom_w,
                       config.dt, -1.0e30f);
  }
  sync_after_launches();
  fast_modes.apply(q1, layout, metrics, config);

  boundary_updater.apply(q1, layout, config.dt);
  compute_slow_tendencies(q1, layout, metrics, config, tendencies);
  for (std::size_t n = 0; n < states.size(); ++n) {
    tmp[n].copy_all_from(q1[n]);
    update_stage_field(tmp[n].rho_d, q1[n].rho_d, tendencies[n].rho_d,
                       config.dt, 1.0e-6f);
    update_stage_field(tmp[n].rho_theta_m, q1[n].rho_theta_m,
                       tendencies[n].rho_theta_m, config.dt, 1.0e-6f);
    update_stage_field(tmp[n].mom_u, q1[n].mom_u, tendencies[n].mom_u,
                       config.dt, -1.0e30f);
    update_stage_field(tmp[n].mom_v, q1[n].mom_v, tendencies[n].mom_v,
                       config.dt, -1.0e30f);
    update_stage_field(tmp[n].mom_w, q1[n].mom_w, tendencies[n].mom_w,
                       config.dt, -1.0e30f);
    blend_stage_field(q2[n].rho_d, q0[n].rho_d, tmp[n].rho_d, 0.75f, 0.25f);
    blend_stage_field(q2[n].rho_theta_m, q0[n].rho_theta_m, tmp[n].rho_theta_m,
                      0.75f, 0.25f);
    blend_stage_field(q2[n].mom_u, q0[n].mom_u, tmp[n].mom_u, 0.75f, 0.25f);
    blend_stage_field(q2[n].mom_v, q0[n].mom_v, tmp[n].mom_v, 0.75f, 0.25f);
    blend_stage_field(q2[n].mom_w, q0[n].mom_w, tmp[n].mom_w, 0.75f, 0.25f);
  }
  sync_after_launches();
  fast_modes.apply(q2, layout, metrics, config);

  boundary_updater.apply(q2, layout, 0.5f * config.dt);
  compute_slow_tendencies(q2, layout, metrics, config, tendencies);
  for (std::size_t n = 0; n < states.size(); ++n) {
    tmp[n].copy_all_from(q2[n]);
    update_stage_field(tmp[n].rho_d, q2[n].rho_d, tendencies[n].rho_d,
                       config.dt, 1.0e-6f);
    update_stage_field(tmp[n].rho_theta_m, q2[n].rho_theta_m,
                       tendencies[n].rho_theta_m, config.dt, 1.0e-6f);
    update_stage_field(tmp[n].mom_u, q2[n].mom_u, tendencies[n].mom_u,
                       config.dt, -1.0e30f);
    update_stage_field(tmp[n].mom_v, q2[n].mom_v, tendencies[n].mom_v,
                       config.dt, -1.0e30f);
    update_stage_field(tmp[n].mom_w, q2[n].mom_w, tendencies[n].mom_w,
                       config.dt, -1.0e30f);
    blend_stage_field(states[n].rho_d, q0[n].rho_d, tmp[n].rho_d,
                      1.0f / 3.0f, 2.0f / 3.0f);
    blend_stage_field(states[n].rho_theta_m, q0[n].rho_theta_m,
                      tmp[n].rho_theta_m, 1.0f / 3.0f, 2.0f / 3.0f);
    blend_stage_field(states[n].mom_u, q0[n].mom_u, tmp[n].mom_u,
                      1.0f / 3.0f, 2.0f / 3.0f);
    blend_stage_field(states[n].mom_v, q0[n].mom_v, tmp[n].mom_v,
                      1.0f / 3.0f, 2.0f / 3.0f);
    blend_stage_field(states[n].mom_w, q0[n].mom_w, tmp[n].mom_w,
                      1.0f / 3.0f, 2.0f / 3.0f);
  }
  sync_after_launches();
  fast_modes.apply(states, layout, metrics, config);
}

}  // namespace gwm::dycore
