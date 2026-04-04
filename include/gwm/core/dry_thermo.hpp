#pragma once

#include <cmath>

#include "gwm/core/types.hpp"

#ifndef __host__
#define __host__
#endif

#ifndef __device__
#define __device__
#endif

namespace gwm::core {

constexpr real kDryGasConstant = 287.05f;
constexpr real kReferencePressure = 100000.0f;
constexpr real kKappa = 0.2854f;

__host__ __device__ inline real dry_theta_from_conserved(real rho_d,
                                                         real rho_theta_m) {
  return rho_theta_m / fmaxf(rho_d, 1.0e-6f);
}

__host__ __device__ inline real dry_pressure_from_rho_theta_m(
    real rho_d, real rho_theta_m) {
  const real theta = dry_theta_from_conserved(rho_d, rho_theta_m);
  const real ratio =
      fmaxf(rho_d * kDryGasConstant * theta / kReferencePressure, 1.0e-6f);
  return kReferencePressure *
         powf(ratio, 1.0f / (1.0f - kKappa));
}

__host__ __device__ inline real dry_rho_from_pressure_theta(real pressure,
                                                            real theta) {
  const real exner =
      powf(fmaxf(pressure / kReferencePressure, 1.0e-6f), kKappa);
  return pressure / fmaxf(kDryGasConstant * theta * exner, 1.0e-6f);
}

__host__ __device__ inline real dry_cp() {
  return kDryGasConstant / kKappa;
}

}  // namespace gwm::core
