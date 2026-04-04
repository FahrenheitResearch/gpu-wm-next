#include "gwm/dycore/dry_transport.hpp"

#include <algorithm>

#include "gwm/comm/halo_exchange.hpp"
#include "gwm/core/cuda_utils.hpp"
#include "gwm/dycore/passive_tracer.hpp"

namespace gwm::dycore {

namespace {

__global__ void upwind_stage_kernel(const real* in, real* out, int nx, int ny,
                                    int nz, int halo, real u_adv, real v_adv,
                                    real dt, real dx, real dy) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  const int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (i >= nx || j >= ny || k >= nz) {
    return;
  }

  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  const auto idx = [&](int ii, int jj, int kk) {
    return (kk * pitch_y + (jj + halo)) * pitch_x + (ii + halo);
  };

  const real q = in[idx(i, j, k)];
  const real dqx =
      (u_adv >= 0.0f) ? (q - in[idx(i - 1, j, k)]) / dx
                      : (in[idx(i + 1, j, k)] - q) / dx;
  const real dqy =
      (v_adv >= 0.0f) ? (q - in[idx(i, j - 1, k)]) / dy
                      : (in[idx(i, j + 1, k)] - q) / dy;

  out[idx(i, j, k)] = q - dt * (u_adv * dqx + v_adv * dqy);
}

__global__ void blend_stage_kernel(real* out, const real* a, const real* b,
                                   int total_size, real wa, real wb) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= total_size) {
    return;
  }
  out[idx] = wa * a[idx] + wb * b[idx];
}

__global__ void passive_tracer_stage_kernel(
    const real* rho_q_in, const real* rho_d, real* rho_q_out,
    const real* mom_u, const real* mom_v, const real* mom_w,
    const real* inv_dz_cell_metric, int metrics_nx, int metrics_ny, int i_begin,
    int j_begin, int nx, int ny, int nz, int halo, real dt, real dx, real dy) {
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
    return (static_cast<std::size_t>(kk) * static_cast<std::size_t>(pitch_y) +
            static_cast<std::size_t>(jj + halo)) *
               static_cast<std::size_t>(nx + 1 + 2 * halo) +
           static_cast<std::size_t>(ii + halo);
  };
  const auto yface_idx = [&](int ii, int jj, int kk) {
    return (static_cast<std::size_t>(kk) *
                static_cast<std::size_t>(ny + 1 + 2 * halo) +
            static_cast<std::size_t>(jj + halo)) *
               static_cast<std::size_t>(pitch_x) +
           static_cast<std::size_t>(ii + halo);
  };
  const auto zface_idx = [&](int ii, int jj, int kk) {
    return (static_cast<std::size_t>(kk) * static_cast<std::size_t>(pitch_y) +
            static_cast<std::size_t>(jj + halo)) *
               static_cast<std::size_t>(pitch_x) +
           static_cast<std::size_t>(ii + halo);
  };
  const auto metric_idx = [&](int ii, int jj, int kk) {
    return (static_cast<std::size_t>(kk) * static_cast<std::size_t>(metrics_ny) +
            static_cast<std::size_t>(jj)) *
               static_cast<std::size_t>(metrics_nx) +
           static_cast<std::size_t>(ii);
  };
  const auto mixing_ratio = [&](int ii, int jj, int kk) {
    const auto index = cell_idx(ii, jj, kk);
    return rho_q_in[index] / fmaxf(rho_d[index], 1.0e-6f);
  };

  const real q_c = mixing_ratio(i, j, k);
  const real q_w = mixing_ratio(i - 1, j, k);
  const real q_e = mixing_ratio(i + 1, j, k);
  const real q_s = mixing_ratio(i, j - 1, k);
  const real q_n = mixing_ratio(i, j + 1, k);
  const real q_b = (k > 0) ? mixing_ratio(i, j, k - 1) : q_c;
  const real q_t = (k + 1 < nz) ? mixing_ratio(i, j, k + 1) : q_c;

  const real mass_flux_w = mom_u[xface_idx(i, j, k)];
  const real mass_flux_e = mom_u[xface_idx(i + 1, j, k)];
  const real mass_flux_s = mom_v[yface_idx(i, j, k)];
  const real mass_flux_n = mom_v[yface_idx(i, j + 1, k)];
  const real mass_flux_b = mom_w[zface_idx(i, j, k)];
  const real mass_flux_t = mom_w[zface_idx(i, j, k + 1)];
  const real inv_dz =
      inv_dz_cell_metric[metric_idx(i_begin + i, j_begin + j, k)];

  const real flux_q_w = mass_flux_w * (mass_flux_w >= 0.0f ? q_w : q_c);
  const real flux_q_e = mass_flux_e * (mass_flux_e >= 0.0f ? q_c : q_e);
  const real flux_q_s = mass_flux_s * (mass_flux_s >= 0.0f ? q_s : q_c);
  const real flux_q_n = mass_flux_n * (mass_flux_n >= 0.0f ? q_c : q_n);
  const real flux_q_b = mass_flux_b * (mass_flux_b >= 0.0f ? q_b : q_c);
  const real flux_q_t = mass_flux_t * (mass_flux_t >= 0.0f ? q_c : q_t);

  const real div_rho_q = (flux_q_e - flux_q_w) / dx +
                         (flux_q_n - flux_q_s) / dy +
                         (flux_q_t - flux_q_b) * inv_dz;
  const auto out_index = cell_idx(i, j, k);
  rho_q_out[out_index] = fmaxf(0.0f, rho_q_in[out_index] - dt * div_rho_q);
}

void run_upwind_stage(const state::Field3D<real>& in, state::Field3D<real>& out,
                      const DryTransportConfig& config) {
  dim3 block(8, 8, 4);
  dim3 grid((in.nx() + block.x - 1) / block.x,
            (in.ny() + block.y - 1) / block.y,
            (in.nz() + block.z - 1) / block.z);

  upwind_stage_kernel<<<grid, block>>>(
      in.data(), out.data(), in.nx(), in.ny(), in.nz(), in.halo(), config.u_adv,
      config.v_adv, config.dt, config.dx, config.dy);
  GWM_CUDA_CHECK(cudaGetLastError());
  GWM_CUDA_CHECK(cudaDeviceSynchronize());
}

void blend_fields(state::Field3D<real>& out, const state::Field3D<real>& a,
                  const state::Field3D<real>& b, real wa, real wb) {
  const int total = static_cast<int>(out.total_size());
  constexpr int block_size = 256;
  const int grid_size = (total + block_size - 1) / block_size;
  blend_stage_kernel<<<grid_size, block_size>>>(out.data(), a.data(), b.data(),
                                                total, wa, wb);
  GWM_CUDA_CHECK(cudaGetLastError());
  GWM_CUDA_CHECK(cudaDeviceSynchronize());
}

std::size_t linear_index_3d(int i, int j, int k, int nx, int ny) {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

void require_tracer_layout_compatibility(
    const std::vector<state::TracerState>& tracers,
    const std::vector<DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout) {
  gwm::require(tracers.size() == dry_states.size(),
               "Tracer/dry-state size mismatch in passive tracer transport");
  gwm::require(tracers.size() == layout.size(),
               "Tracer/layout size mismatch in passive tracer transport");
  if (tracers.empty()) {
    return;
  }

  const auto tracer_count = tracers.front().size();
  for (std::size_t n = 0; n < tracers.size(); ++n) {
    gwm::require(tracers[n].size() == tracer_count,
                 "Tracer count must match across layout");
    gwm::require(tracers[n].nx() == dry_states[n].rho_d.nx() &&
                     tracers[n].ny() == dry_states[n].rho_d.ny() &&
                     tracers[n].nz() == dry_states[n].rho_d.nz() &&
                     tracers[n].halo() == dry_states[n].rho_d.halo(),
                 "Tracer field shape mismatch against dry state");
    for (std::size_t tracer = 0; tracer < tracer_count; ++tracer) {
      gwm::require(
          tracers[n].registry().at(static_cast<int>(tracer)).name ==
              tracers.front().registry().at(static_cast<int>(tracer)).name,
          "Tracer registry must match across layout");
    }
  }
}

void prepare_dry_transport_halos(
    const std::vector<DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    std::vector<state::Field3D<real>>& rho_fields,
    std::vector<state::FaceField<real>>& mom_u_fields,
    std::vector<state::FaceField<real>>& mom_v_fields,
    std::vector<state::FaceField<real>>& mom_w_fields) {
  rho_fields.clear();
  mom_u_fields.clear();
  mom_v_fields.clear();
  mom_w_fields.clear();
  rho_fields.reserve(dry_states.size());
  mom_u_fields.reserve(dry_states.size());
  mom_v_fields.reserve(dry_states.size());
  mom_w_fields.reserve(dry_states.size());

  for (const auto& state : dry_states) {
    rho_fields.push_back(state.rho_d.clone_empty_like("_rho_tracer_halo"));
    rho_fields.back().copy_all_from(state.rho_d);
    mom_u_fields.emplace_back(state.rho_d.nx(), state.rho_d.ny(),
                              state.rho_d.nz(), state.rho_d.halo(),
                              state::FaceOrientation::X, "mom_u_tracer_halo");
    mom_u_fields.back().storage().copy_all_from(state.mom_u.storage());
    mom_v_fields.emplace_back(state.rho_d.nx(), state.rho_d.ny(),
                              state.rho_d.nz(), state.rho_d.halo(),
                              state::FaceOrientation::Y, "mom_v_tracer_halo");
    mom_v_fields.back().storage().copy_all_from(state.mom_v.storage());
    mom_w_fields.emplace_back(state.rho_d.nx(), state.rho_d.ny(),
                              state.rho_d.nz(), state.rho_d.halo(),
                              state::FaceOrientation::Z, "mom_w_tracer_halo");
    mom_w_fields.back().storage().copy_all_from(state.mom_w.storage());
  }

  comm::HaloExchange::exchange_scalar(rho_fields, layout);
  comm::HaloExchange::exchange_face(mom_u_fields, layout);
  comm::HaloExchange::exchange_face(mom_v_fields, layout);
  comm::HaloExchange::exchange_face(mom_w_fields, layout);
  comm::HaloExchange::synchronize_owned_face_interfaces(mom_u_fields, layout);
  comm::HaloExchange::synchronize_owned_face_interfaces(mom_v_fields, layout);
  comm::HaloExchange::exchange_face(mom_u_fields, layout);
  comm::HaloExchange::exchange_face(mom_v_fields, layout);
}

void run_passive_tracer_stage(
    const state::Field3D<real>& rho_q_in, const state::Field3D<real>& rho_d,
    const state::FaceField<real>& mom_u, const state::FaceField<real>& mom_v,
    const state::FaceField<real>& mom_w, state::Field3D<real>& rho_q_out,
    const domain::SubdomainDescriptor& desc, const domain::GridMetrics& metrics,
    real dt) {
  dim3 block(8, 8, 4);
  dim3 grid((rho_q_in.nx() + block.x - 1) / block.x,
            (rho_q_in.ny() + block.y - 1) / block.y,
            (rho_q_in.nz() + block.z - 1) / block.z);

  passive_tracer_stage_kernel<<<grid, block>>>(
      rho_q_in.data(), rho_d.data(), rho_q_out.data(), mom_u.storage().data(),
      mom_v.storage().data(), mom_w.storage().data(),
      metrics.inv_dz_cell_metric.data(), metrics.nx, metrics.ny, desc.i_begin,
      desc.j_begin, rho_q_in.nx(), rho_q_in.ny(), rho_q_in.nz(),
      rho_q_in.halo(), dt, static_cast<real>(metrics.dx),
      static_cast<real>(metrics.dy));
  GWM_CUDA_CHECK(cudaGetLastError());
  GWM_CUDA_CHECK(cudaDeviceSynchronize());
}

void apply_tracer_boundaries_to_fields(
    std::vector<state::Field3D<real>>& fields,
    const std::vector<state::TracerState>& tracer_layout,
    std::size_t tracer_index,
    const std::vector<domain::SubdomainDescriptor>& layout, real sim_time,
    TracerBoundaryUpdater& boundary_updater) {
  std::vector<state::TracerState> staged_tracers;
  staged_tracers.reserve(fields.size());
  const auto tracer_spec =
      tracer_layout.front().registry().at(static_cast<int>(tracer_index));
  for (std::size_t n = 0; n < fields.size(); ++n) {
    state::TracerRegistry registry;
    registry.add(tracer_spec);
    staged_tracers.emplace_back(std::move(registry), fields[n].nx(), fields[n].ny(),
                                fields[n].nz(), fields[n].halo(),
                                "stage_tracer");
    staged_tracers.back().mass(0).copy_all_from(fields[n]);
  }
  boundary_updater.apply(staged_tracers, layout, sim_time);
  for (std::size_t n = 0; n < fields.size(); ++n) {
    fields[n].copy_all_from(staged_tracers[n].mass(0));
  }
}

}  // namespace

void NullTracerBoundaryUpdater::apply(
    std::vector<state::TracerState>& /*states*/,
    const std::vector<domain::SubdomainDescriptor>& /*layout*/,
    real /*sim_time*/) {}

void advance_scalar_ssprk3(
    std::vector<state::Field3D<real>>& fields,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const DryTransportConfig& config) {
  gwm::require(fields.size() == layout.size(),
               "Field/layout size mismatch in advance_scalar_ssprk3");

  std::vector<state::Field3D<real>> q0;
  std::vector<state::Field3D<real>> q1;
  std::vector<state::Field3D<real>> q2;
  std::vector<state::Field3D<real>> tmp;
  q0.reserve(fields.size());
  q1.reserve(fields.size());
  q2.reserve(fields.size());
  tmp.reserve(fields.size());

  for (const auto& field : fields) {
    q0.push_back(field.clone_empty_like("_q0"));
    q0.back().copy_all_from(field);
    q1.push_back(field.clone_empty_like("_q1"));
    q2.push_back(field.clone_empty_like("_q2"));
    tmp.push_back(field.clone_empty_like("_tmp"));
  }

  comm::HaloExchange::exchange_scalar(fields, layout);
  for (std::size_t n = 0; n < fields.size(); ++n) {
    run_upwind_stage(fields[n], q1[n], config);
  }

  comm::HaloExchange::exchange_scalar(q1, layout);
  for (std::size_t n = 0; n < fields.size(); ++n) {
    run_upwind_stage(q1[n], tmp[n], config);
    blend_fields(q2[n], q0[n], tmp[n], 0.75f, 0.25f);
  }

  comm::HaloExchange::exchange_scalar(q2, layout);
  for (std::size_t n = 0; n < fields.size(); ++n) {
    run_upwind_stage(q2[n], tmp[n], config);
    blend_fields(fields[n], q0[n], tmp[n], 1.0f / 3.0f, 2.0f / 3.0f);
  }
}

std::vector<state::TracerState> make_specific_humidity_tracers_from_global_field(
    const std::vector<DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout, int global_nx,
    int global_ny, const std::vector<real>& specific_humidity,
    const std::string& label_prefix) {
  gwm::require(dry_states.size() == layout.size(),
               "Dry-state/layout size mismatch in "
               "make_specific_humidity_tracers_from_global_field");
  gwm::require(!layout.empty(),
               "Layout must not be empty in "
               "make_specific_humidity_tracers_from_global_field");
  const auto global_nz = layout.front().nz;
  gwm::require(specific_humidity.size() ==
                   static_cast<std::size_t>(global_nx) *
                       static_cast<std::size_t>(global_ny) *
                       static_cast<std::size_t>(global_nz),
               "Specific-humidity field size mismatch");

  std::vector<state::TracerState> tracers;
  tracers.reserve(layout.size());
  for (const auto& desc : layout) {
    state::TracerState tracer_state(state::make_specific_humidity_registry(),
                                    desc.nx_local(), desc.ny_local(), desc.nz,
                                    desc.halo,
                                    label_prefix + "_rank_" +
                                        std::to_string(desc.rank));
    tracer_state.fill_zero();
    auto& rho_q = tracer_state.mass(state::kSpecificHumidityTracerName);
    for (int j = 0; j < desc.ny_local(); ++j) {
      const int j_global = desc.j_begin + j;
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k = 0; k < desc.nz; ++k) {
          const auto global_index =
              linear_index_3d(i_global, j_global, k, global_nx, global_ny);
          rho_q(i, j, k) =
              dry_states[static_cast<std::size_t>(desc.rank)].rho_d(i, j, k) *
              specific_humidity[global_index];
        }
      }
    }
    tracers.push_back(std::move(tracer_state));
  }
  return tracers;
}

std::vector<state::TracerState> make_warm_rain_tracers_from_global_fields(
    const std::vector<DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout, int global_nx,
    int global_ny, const std::vector<real>& specific_humidity,
    const std::vector<real>* cloud_water, const std::vector<real>* rain_water,
    const std::string& label_prefix) {
  gwm::require(dry_states.size() == layout.size(),
               "Dry-state/layout size mismatch in "
               "make_warm_rain_tracers_from_global_fields");
  gwm::require(!layout.empty(),
               "Layout must not be empty in "
               "make_warm_rain_tracers_from_global_fields");
  const auto global_nz = layout.front().nz;
  const auto expected_size = static_cast<std::size_t>(global_nx) *
                             static_cast<std::size_t>(global_ny) *
                             static_cast<std::size_t>(global_nz);
  gwm::require(specific_humidity.size() == expected_size,
               "Specific-humidity field size mismatch");
  if (cloud_water != nullptr) {
    gwm::require(cloud_water->size() == expected_size,
                 "Cloud-water field size mismatch");
  }
  if (rain_water != nullptr) {
    gwm::require(rain_water->size() == expected_size,
                 "Rain-water field size mismatch");
  }

  std::vector<state::TracerState> tracers;
  tracers.reserve(layout.size());
  for (const auto& desc : layout) {
    state::TracerState tracer_state(state::make_warm_rain_registry(),
                                    desc.nx_local(), desc.ny_local(), desc.nz,
                                    desc.halo,
                                    label_prefix + "_rank_" +
                                        std::to_string(desc.rank));
    tracer_state.fill_zero();
    auto& rho_qv = tracer_state.mass(state::kSpecificHumidityTracerName);
    auto& rho_qc = tracer_state.mass(state::kCloudWaterTracerName);
    auto& rho_qr = tracer_state.mass(state::kRainWaterTracerName);
    const auto& rho_d =
        dry_states[static_cast<std::size_t>(desc.rank)].rho_d;
    for (int j = 0; j < desc.ny_local(); ++j) {
      const int j_global = desc.j_begin + j;
      for (int i = 0; i < desc.nx_local(); ++i) {
        const int i_global = desc.i_begin + i;
        for (int k = 0; k < desc.nz; ++k) {
          const auto global_index =
              linear_index_3d(i_global, j_global, k, global_nx, global_ny);
          const real rho = rho_d(i, j, k);
          rho_qv(i, j, k) =
              rho * std::max<real>(0.0f, specific_humidity[global_index]);
          rho_qc(i, j, k) =
              rho * std::max<real>(0.0f, cloud_water == nullptr
                                              ? 0.0f
                                              : (*cloud_water)[global_index]);
          rho_qr(i, j, k) =
              rho * std::max<real>(0.0f, rain_water == nullptr
                                              ? 0.0f
                                              : (*rain_water)[global_index]);
        }
      }
    }
    tracers.push_back(std::move(tracer_state));
  }
  return tracers;
}

void advance_passive_tracers_ssprk3(
    std::vector<state::TracerState>& tracers,
    const std::vector<DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const DryStepperConfig& config,
    TracerBoundaryUpdater& boundary_updater) {
  require_tracer_layout_compatibility(tracers, dry_states, layout);
  if (tracers.empty() || tracers.front().size() == 0) {
    return;
  }

  std::vector<state::Field3D<real>> rho_fields;
  std::vector<state::FaceField<real>> mom_u_fields;
  std::vector<state::FaceField<real>> mom_v_fields;
  std::vector<state::FaceField<real>> mom_w_fields;
  prepare_dry_transport_halos(dry_states, layout, rho_fields, mom_u_fields,
                              mom_v_fields, mom_w_fields);

  for (std::size_t tracer_index = 0; tracer_index < tracers.front().size();
       ++tracer_index) {
    std::vector<state::Field3D<real>> q0;
    std::vector<state::Field3D<real>> q1;
    std::vector<state::Field3D<real>> q2;
    std::vector<state::Field3D<real>> tmp;
    q0.reserve(tracers.size());
    q1.reserve(tracers.size());
    q2.reserve(tracers.size());
    tmp.reserve(tracers.size());

    for (const auto& tracer_state : tracers) {
      const auto& field = tracer_state.mass(static_cast<int>(tracer_index));
      q0.push_back(field.clone_empty_like("_q0"));
      q0.back().copy_all_from(field);
      q1.push_back(field.clone_empty_like("_q1"));
      q2.push_back(field.clone_empty_like("_q2"));
      tmp.push_back(field.clone_empty_like("_tmp"));
    }

    std::vector<state::Field3D<real>> stage_fields;
    stage_fields.reserve(tracers.size());
    for (const auto& tracer_state : tracers) {
      stage_fields.push_back(
          tracer_state.mass(static_cast<int>(tracer_index)).clone_empty_like(
              "_stage"));
      stage_fields.back().copy_all_from(
          tracer_state.mass(static_cast<int>(tracer_index)));
    }

    apply_tracer_boundaries_to_fields(stage_fields, tracers, tracer_index, layout,
                                      0.0f, boundary_updater);
    comm::HaloExchange::exchange_scalar(stage_fields, layout);
    for (std::size_t n = 0; n < tracers.size(); ++n) {
      run_passive_tracer_stage(stage_fields[n], rho_fields[n], mom_u_fields[n],
                               mom_v_fields[n], mom_w_fields[n], q1[n],
                               layout[n], metrics, config.dt);
    }

    apply_tracer_boundaries_to_fields(q1, tracers, tracer_index, layout,
                                      config.dt, boundary_updater);
    comm::HaloExchange::exchange_scalar(q1, layout);
    for (std::size_t n = 0; n < tracers.size(); ++n) {
      run_passive_tracer_stage(q1[n], rho_fields[n], mom_u_fields[n],
                               mom_v_fields[n], mom_w_fields[n], tmp[n],
                               layout[n], metrics, config.dt);
      blend_fields(q2[n], q0[n], tmp[n], 0.75f, 0.25f);
    }

    apply_tracer_boundaries_to_fields(q2, tracers, tracer_index, layout,
                                      0.5f * config.dt, boundary_updater);
    comm::HaloExchange::exchange_scalar(q2, layout);
    for (std::size_t n = 0; n < tracers.size(); ++n) {
      run_passive_tracer_stage(q2[n], rho_fields[n], mom_u_fields[n],
                               mom_v_fields[n], mom_w_fields[n], tmp[n],
                               layout[n], metrics, config.dt);
      blend_fields(tracers[n].mass(static_cast<int>(tracer_index)), q0[n],
                   tmp[n], 1.0f / 3.0f, 2.0f / 3.0f);
    }
  }
}

}  // namespace gwm::dycore
