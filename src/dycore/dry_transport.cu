#include "gwm/dycore/dry_transport.hpp"

#include "gwm/comm/halo_exchange.hpp"
#include "gwm/core/cuda_utils.hpp"

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

}  // namespace

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

}  // namespace gwm::dycore
