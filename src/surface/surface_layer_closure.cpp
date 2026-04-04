#include "gwm/surface/surface_layer_closure.hpp"

#include <algorithm>
#include <cmath>

#include "gwm/core/moist_thermo.hpp"

namespace gwm::surface {

namespace {

constexpr real kMinRoughness = 1.0e-4f;
constexpr real kMinPressurePa = 1000.0f;
constexpr real kVonKarman = 0.4f;

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
  return core::saturation_specific_humidity(temperature_k, pressure_pa);
}

real relative_humidity_from_specific_humidity(real specific_humidity,
                                              real temperature_k,
                                              real pressure_pa) {
  return core::relative_humidity_from_specific_humidity(
      specific_humidity, temperature_k, pressure_pa);
}

real dewpoint_from_specific_humidity(real specific_humidity,
                                     real pressure_pa) {
  return core::dewpoint_from_specific_humidity(specific_humidity, pressure_pa);
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
