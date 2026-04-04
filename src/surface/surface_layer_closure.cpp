#include "gwm/surface/surface_layer_closure.hpp"

#include <algorithm>
#include <cmath>

namespace gwm::surface {

namespace {

constexpr real kMinRoughness = 1.0e-4f;
constexpr real kMinPressurePa = 1000.0f;
constexpr real kVonKarman = 0.4f;
constexpr real kRvOverRd = 0.622f;

real clamp_height(real z_target, real z0) {
  return std::max(z_target, z0 * 1.01f);
}

real interpolation_weight(real z_target, real z_ref, real z0) {
  const real z_target_safe = clamp_height(z_target, z0);
  const real z_ref_safe = clamp_height(z_ref, z0);
  const real denom = std::log(z_ref_safe / z0);
  return (denom > 0.0f) ? (std::log(z_target_safe / z0) / denom) : 1.0f;
}

}  // namespace

real saturation_specific_humidity(real temperature_k, real pressure_pa) {
  const real temp_c = temperature_k - 273.15f;
  const real es_hpa =
      6.112f * std::exp((17.67f * temp_c) / (temp_c + 243.5f));
  const real es_pa = es_hpa * 100.0f;
  const real denom = std::max(pressure_pa - es_pa, kMinPressurePa);
  return std::clamp(kRvOverRd * es_pa / denom, 0.0f, 0.05f);
}

real relative_humidity_from_specific_humidity(real specific_humidity,
                                              real temperature_k,
                                              real pressure_pa) {
  const real qsat = saturation_specific_humidity(temperature_k, pressure_pa);
  return std::clamp(
      (qsat > 0.0f) ? (100.0f * specific_humidity / qsat) : 0.0f, 0.0f, 100.0f);
}

real dewpoint_from_specific_humidity(real specific_humidity,
                                     real pressure_pa) {
  const real q = std::clamp(specific_humidity, 0.0f, 1.0f);
  const real vapor_pressure_pa =
      pressure_pa * q / std::max(kRvOverRd + (1.0f - kRvOverRd) * q, 1.0e-6f);
  const real vapor_pressure_hpa =
      std::max(vapor_pressure_pa * 0.01f, 1.0e-6f);
  const real log_term = std::log(vapor_pressure_hpa / 6.112f);
  const real dewpoint_c = (243.5f * log_term) / (17.67f - log_term);
  return dewpoint_c + 273.15f;
}

SurfaceLayerDiagnostics evaluate_neutral_surface_layer(
    const SurfaceLayerInputs& in) {
  gwm::require(in.z_ref > 0.0f, "Surface-layer closure requires z_ref > 0");
  gwm::require(in.z_target_temp > 0.0f,
               "Surface-layer closure requires z_target_temp > 0");
  gwm::require(in.z_target_wind > 0.0f,
               "Surface-layer closure requires z_target_wind > 0");
  gwm::require(in.psfc > kMinPressurePa,
               "Surface-layer closure requires physical psfc");

  const real z0m = std::max(in.z0m, kMinRoughness);
  const real z0h = std::max(in.z0h, kMinRoughness);
  const real wind_weight = interpolation_weight(in.z_target_wind, in.z_ref, z0m);
  const real scalar_weight = interpolation_weight(in.z_target_temp, in.z_ref, z0h);

  SurfaceLayerDiagnostics out{};
  out.u10 = in.u_ref * wind_weight;
  out.v10 = in.v_ref * wind_weight;
  out.t2 = in.tskin + (in.theta_ref - in.tskin) * scalar_weight;
  out.q2 = in.q_surface + (in.q_ref - in.q_surface) * scalar_weight;
  out.q2 = std::clamp(out.q2, 0.0f, std::max(in.q_ref, in.q_surface));

  out.wind_speed_ref =
      std::sqrt(in.u_ref * in.u_ref + in.v_ref * in.v_ref);
  out.wind_speed_target = std::sqrt(out.u10 * out.u10 + out.v10 * out.v10);

  const real log_m = std::log(clamp_height(in.z_ref, z0m) / z0m);
  const real log_h = std::log(clamp_height(in.z_ref, z0h) / z0h);
  out.cm = (log_m > 0.0f) ? (kVonKarman * kVonKarman) / (log_m * log_m) : 0.0f;
  out.ch = (log_m > 0.0f && log_h > 0.0f)
               ? (kVonKarman * kVonKarman) / (log_m * log_h)
               : 0.0f;
  out.cq = out.ch;
  out.ustar = std::sqrt(std::max(out.cm, 0.0f)) * out.wind_speed_ref;

  const real qsat_2 = saturation_specific_humidity(out.t2, in.psfc);
  out.rh2 = std::clamp((qsat_2 > 0.0f) ? (100.0f * out.q2 / qsat_2) : 0.0f,
                       0.0f, 100.0f);
  return out;
}

}  // namespace gwm::surface
