#include "gwm/physics/warm_rain.hpp"

#include <algorithm>
#include <sstream>

#include "gwm/core/cuda_utils.hpp"
#include "gwm/core/moist_thermo.hpp"
#include "gwm/state/tracer_registry.hpp"

namespace gwm::physics {

namespace {

using gwm::state::kCloudWaterTracerName;
using gwm::state::kRainWaterTracerName;
using gwm::state::kSpecificHumidityTracerName;

__global__ void warm_rain_kernel(
    const real* rho_d, real* rho_theta_m, real* rho_qv, real* rho_qc,
    real* rho_qr, int nx, int ny, int nz, int halo, real dt,
    real condensation_relaxation, real cloud_autoconversion_threshold,
    real cloud_autoconversion_rate, real rain_evaporation_rate,
    bool enable_latent_heating) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  const int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (i >= nx || j >= ny || k >= nz) {
    return;
  }

  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  const auto idx = (k * pitch_y + (j + halo)) * pitch_x + (i + halo);

  const real rho = fmaxf(rho_d[idx], 1.0e-6f);
  real theta_m = rho_theta_m[idx] / rho;
  real pressure = core::dry_pressure_from_rho_theta_m(rho, rho_theta_m[idx]);
  real temperature = core::dry_temperature_from_conserved(rho, rho_theta_m[idx]);
  const real exner =
      fmaxf(temperature / fmaxf(theta_m, 1.0e-6f), 1.0e-6f);

  real qv = fmaxf(rho_qv[idx] / rho, 0.0f);
  real qc = fmaxf(rho_qc[idx] / rho, 0.0f);
  real qr = fmaxf(rho_qr[idx] / rho, 0.0f);

  real qsat = core::saturation_specific_humidity(temperature, pressure);
  if (qv > qsat) {
    const real condense =
        fminf((qv - qsat) * fmaxf(condensation_relaxation, 0.0f), qv);
    qv -= condense;
    qc += condense;
    if (enable_latent_heating && condense > 0.0f) {
      theta_m += (core::kLatentHeatVaporization /
                  (core::dry_cp() * fmaxf(exner, 1.0e-6f))) *
                 condense;
    }
  } else if (qv < qsat && qc > 0.0f) {
    const real evaporate =
        fminf((qsat - qv) * fmaxf(condensation_relaxation, 0.0f), qc);
    qv += evaporate;
    qc -= evaporate;
    if (enable_latent_heating && evaporate > 0.0f) {
      theta_m -= (core::kLatentHeatVaporization /
                  (core::dry_cp() * fmaxf(exner, 1.0e-6f))) *
                 evaporate;
    }
  }

  const real excess_cloud = fmaxf(qc - cloud_autoconversion_threshold, 0.0f);
  if (excess_cloud > 0.0f && cloud_autoconversion_rate > 0.0f) {
    const real convert =
        fminf(excess_cloud, dt * cloud_autoconversion_rate * excess_cloud);
    qc -= convert;
    qr += convert;
  }

  pressure = core::dry_pressure_from_rho_theta_m(rho, rho * theta_m);
  temperature = theta_m *
                powf(fmaxf(pressure / core::kReferencePressure, 1.0e-6f),
                     core::kKappa);
  qsat = core::saturation_specific_humidity(temperature, pressure);
  if (qr > 0.0f && qv < qsat && rain_evaporation_rate > 0.0f) {
    const real evaporate =
        fminf(qr, dt * rain_evaporation_rate * fminf(qsat - qv, qr));
    qr -= evaporate;
    qv += evaporate;
    if (enable_latent_heating && evaporate > 0.0f) {
      theta_m -= (core::kLatentHeatVaporization /
                  (core::dry_cp() * fmaxf(exner, 1.0e-6f))) *
                 evaporate;
    }
  }

  qv = fmaxf(qv, 0.0f);
  qc = fmaxf(qc, 0.0f);
  qr = fmaxf(qr, 0.0f);

  rho_qv[idx] = rho * qv;
  rho_qc[idx] = rho * qc;
  rho_qr[idx] = rho * qr;
  rho_theta_m[idx] = rho * theta_m;
}

void require_warm_rain_tracer_contract(const std::vector<dycore::DryState>& states,
                                       const std::vector<state::TracerState>& tracers) {
  gwm::require(states.size() == tracers.size(),
               "State/tracer size mismatch in warm-rain microphysics");
  for (std::size_t n = 0; n < tracers.size(); ++n) {
    gwm::require(tracers[n].find(kSpecificHumidityTracerName).has_value(),
                 "Warm-rain microphysics requires specific_humidity tracer");
    gwm::require(tracers[n].find(kCloudWaterTracerName).has_value(),
                 "Warm-rain microphysics requires cloud_water_mixing_ratio tracer");
    gwm::require(tracers[n].find(kRainWaterTracerName).has_value(),
                 "Warm-rain microphysics requires rain_water_mixing_ratio tracer");
    gwm::require(tracers[n].nx() == states[n].rho_d.nx() &&
                     tracers[n].ny() == states[n].rho_d.ny() &&
                     tracers[n].nz() == states[n].rho_d.nz() &&
                     tracers[n].halo() == states[n].rho_d.halo(),
                 "Warm-rain tracer fields must match dry-state shape");
  }
}

}  // namespace

void apply_warm_rain_microphysics(std::vector<dycore::DryState>& states,
                                  std::vector<state::TracerState>& tracers,
                                  const WarmRainConfig& config) {
  require_warm_rain_tracer_contract(states, tracers);
  if (states.empty()) {
    return;
  }

  dim3 block(8, 8, 4);
  for (std::size_t n = 0; n < states.size(); ++n) {
    auto& state = states[n];
    auto& tracer = tracers[n];
    const auto qv_index = *tracer.find(kSpecificHumidityTracerName);
    const auto qc_index = *tracer.find(kCloudWaterTracerName);
    const auto qr_index = *tracer.find(kRainWaterTracerName);

    dim3 grid((state.rho_d.nx() + block.x - 1) / block.x,
              (state.rho_d.ny() + block.y - 1) / block.y,
              (state.rho_d.nz() + block.z - 1) / block.z);

    warm_rain_kernel<<<grid, block>>>(
        state.rho_d.data(), state.rho_theta_m.data(),
        tracer.mass(qv_index).data(), tracer.mass(qc_index).data(),
        tracer.mass(qr_index).data(), state.rho_d.nx(), state.rho_d.ny(),
        state.rho_d.nz(), state.rho_d.halo(), config.dt,
        config.condensation_relaxation,
        config.cloud_autoconversion_threshold, config.cloud_autoconversion_rate,
        config.rain_evaporation_rate, config.enable_latent_heating);
    GWM_CUDA_CHECK(cudaGetLastError());
  }
  GWM_CUDA_CHECK(cudaDeviceSynchronize());
}

WarmRainSummary summarize_warm_rain(
    const std::vector<state::TracerState>& tracers) {
  WarmRainSummary summary{};
  for (const auto& tracer : tracers) {
    if (const auto qv_index = tracer.find(kSpecificHumidityTracerName);
        qv_index.has_value()) {
      summary.total_water_vapor_mass += tracer.total_mass(*qv_index);
    }
    if (const auto qc_index = tracer.find(kCloudWaterTracerName);
        qc_index.has_value()) {
      summary.total_cloud_water_mass += tracer.total_mass(*qc_index);
    }
    if (const auto qr_index = tracer.find(kRainWaterTracerName);
        qr_index.has_value()) {
      summary.total_rain_water_mass += tracer.total_mass(*qr_index);
    }
  }
  return summary;
}

std::string warm_rain_summary_to_json(const WarmRainSummary& summary,
                                      const std::string& indent) {
  std::ostringstream oss;
  oss << indent << "{\n";
  oss << indent << "  \"total_water_vapor_mass\": "
      << summary.total_water_vapor_mass << ",\n";
  oss << indent << "  \"total_cloud_water_mass\": "
      << summary.total_cloud_water_mass << ",\n";
  oss << indent << "  \"total_rain_water_mass\": "
      << summary.total_rain_water_mass << "\n";
  oss << indent << "}";
  return oss.str();
}

}  // namespace gwm::physics
