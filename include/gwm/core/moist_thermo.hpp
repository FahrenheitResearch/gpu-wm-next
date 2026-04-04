#pragma once

#include <cmath>

#include "gwm/core/dry_thermo.hpp"
#include "gwm/core/types.hpp"

#ifndef __host__
#define __host__
#endif

#ifndef __device__
#define __device__
#endif

namespace gwm::core {

constexpr real kWaterVaporGasConstant = 461.5f;
constexpr real kLatentHeatVaporization = 2.5e6f;
constexpr real kRvOverRd = kWaterVaporGasConstant / kDryGasConstant;
constexpr real kEpsilon = kDryGasConstant / kWaterVaporGasConstant;

__host__ __device__ inline real saturation_vapor_pressure_pa(
    real temperature_k) {
  const real temp_c = temperature_k - 273.15f;
  const real es_hpa =
      6.112f * expf((17.67f * temp_c) / fmaxf(temp_c + 243.5f, 1.0e-6f));
  return es_hpa * 100.0f;
}

__host__ __device__ inline real saturation_specific_humidity(
    real temperature_k, real pressure_pa) {
  const real es_pa = saturation_vapor_pressure_pa(temperature_k);
  const real denom = fmaxf(pressure_pa - es_pa, 1000.0f);
  return fminf(fmaxf(kEpsilon * es_pa / denom, 0.0f), 0.05f);
}

__host__ __device__ inline real relative_humidity_from_specific_humidity(
    real specific_humidity, real temperature_k, real pressure_pa) {
  const real qsat = saturation_specific_humidity(temperature_k, pressure_pa);
  return fminf(fmaxf((qsat > 0.0f) ? (100.0f * specific_humidity / qsat) : 0.0f,
                     0.0f),
               100.0f);
}

__host__ __device__ inline real dewpoint_from_specific_humidity(
    real specific_humidity, real pressure_pa) {
  const real q = fminf(fmaxf(specific_humidity, 0.0f), 1.0f);
  const real vapor_pressure_pa =
      pressure_pa * q / fmaxf(kEpsilon + (1.0f - kEpsilon) * q, 1.0e-6f);
  const real vapor_pressure_hpa = fmaxf(vapor_pressure_pa * 0.01f, 1.0e-6f);
  const real log_term = logf(vapor_pressure_hpa / 6.112f);
  const real dewpoint_c = (243.5f * log_term) / (17.67f - log_term);
  return dewpoint_c + 273.15f;
}

__host__ __device__ inline real dry_temperature_from_conserved(real rho_d,
                                                               real rho_theta_m) {
  const real theta = dry_theta_from_conserved(rho_d, rho_theta_m);
  const real pressure = dry_pressure_from_rho_theta_m(rho_d, rho_theta_m);
  const real exner =
      powf(fmaxf(pressure / kReferencePressure, 1.0e-6f), kKappa);
  return theta * exner;
}

__host__ __device__ inline real synthetic_reflectivity_dbz_from_rain(
    real rho_d, real rain_mixing_ratio) {
  const real qr = fmaxf(rain_mixing_ratio, 0.0f);
  if (qr <= 1.0e-8f || rho_d <= 1.0e-8f) {
    return -20.0f;
  }
  const real z_linear = 3.63e9f * powf(rho_d * qr, 1.75f);
  return 10.0f * log10f(fmaxf(z_linear, 1.0e-6f));
}

}  // namespace gwm::core
