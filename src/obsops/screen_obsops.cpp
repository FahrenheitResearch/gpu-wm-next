#include "gwm/obsops/screen_obsops.hpp"

namespace gwm::obsops {

namespace {

constexpr real kMinRoughness = 1.0e-4f;
constexpr real kRvOverRd = 0.622f;

real log_interp_weight(real z_target, real z_ref, real z0) {
  const real zr = std::max(z_ref, z0 * 1.01f);
  const real zt = std::max(z_target, z0 * 1.01f);
  const real denom = std::log(zr / z0);
  return (denom > 0.0f) ? (std::log(zt / z0) / denom) : 1.0f;
}

}  // namespace

real saturation_specific_humidity(real temperature_k, real pressure_pa) {
  const real temp_c = temperature_k - 273.15f;
  const real es_hpa =
      6.112f * std::exp((17.67f * temp_c) / (temp_c + 243.5f));
  const real es_pa = es_hpa * 100.0f;
  return std::clamp(kRvOverRd * es_pa / std::max(pressure_pa - es_pa, 1000.0f),
                    0.0f, 0.05f);
}

ScreenDiagnostics diagnose_neutral(const ScreenInputs& in) {
  const real z0m = std::max(in.z0m, kMinRoughness);
  const real z0h = std::max(in.z0h, kMinRoughness);
  const real wind_weight = log_interp_weight(in.z_target_wind, in.z_ref, z0m);
  const real scalar_weight =
      log_interp_weight(in.z_target_temp, in.z_ref, z0h);

  ScreenDiagnostics out{};
  out.u10 = in.u_ref * wind_weight;
  out.v10 = in.v_ref * wind_weight;
  out.t2 = in.tskin + (in.theta_ref - in.tskin) * scalar_weight;
  out.q2 = in.q_surface + (in.q_ref - in.q_surface) * scalar_weight;
  out.q2 = std::max(out.q2, 0.0f);

  const real qsat_2 = saturation_specific_humidity(out.t2, in.psfc);
  out.rh2 = std::clamp((qsat_2 > 0.0f) ? 100.0f * out.q2 / qsat_2 : 0.0f, 0.0f,
                       100.0f);
  return out;
}

}  // namespace gwm::obsops
