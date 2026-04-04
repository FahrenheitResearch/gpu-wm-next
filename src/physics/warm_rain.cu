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

__global__ void rain_fallout_nominal_kernel(
    real* rho_qr, real* accumulated_surface_precip, int nx, int ny, int nz,
    int halo, real dt, real terminal_velocity, real layer_depth_m) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  if (i >= nx || j >= ny) {
    return;
  }

  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  const real dz = fmaxf(layer_depth_m, 1.0e-3f);
  const real fallout_fraction =
      fminf(fmaxf(terminal_velocity * dt / dz, 0.0f), 1.0f);
  if (fallout_fraction <= 0.0f) {
    return;
  }

  real incoming_mass_per_area = 0.0f;
  for (int k = nz - 1; k >= 0; --k) {
    const auto idx = (k * pitch_y + (j + halo)) * pitch_x + (i + halo);
    const real old_mass_per_volume = fmaxf(rho_qr[idx], 0.0f);
    const real old_mass_per_area = old_mass_per_volume * dz;
    const real outgoing_mass_per_area = fallout_fraction * old_mass_per_area;
    const real retained_mass_per_area = old_mass_per_area - outgoing_mass_per_area;
    const real updated_mass_per_area =
        retained_mass_per_area + incoming_mass_per_area;
    rho_qr[idx] = updated_mass_per_area / dz;
    incoming_mass_per_area = outgoing_mass_per_area;
  }
  accumulated_surface_precip[static_cast<std::size_t>(j) * nx + i] +=
      incoming_mass_per_area;
}

__global__ void rain_fallout_metric_kernel(
    const real* rho_qr_old, real* rho_qr, real* accumulated_surface_precip,
    int nx, int ny, int nz, int halo, const real* inv_dz_metric,
    int metrics_nx, int metrics_ny, int i_begin, int j_begin, real dt,
    real terminal_velocity) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;
  if (i >= nx || j >= ny) {
    return;
  }

  const int pitch_x = nx + 2 * halo;
  const int pitch_y = ny + 2 * halo;
  real incoming_mass_per_area = 0.0f;

  for (int k = nz - 1; k >= 0; --k) {
    const auto local_idx = (k * pitch_y + (j + halo)) * pitch_x + (i + halo);
    const auto metric_idx =
        (k * metrics_ny + (j_begin + j)) * metrics_nx + (i_begin + i);
    const real inv_dz = fmaxf(inv_dz_metric[metric_idx], 1.0e-6f);
    const real dz = 1.0f / inv_dz;
    const real fallout_fraction =
        fminf(fmaxf(terminal_velocity * dt * inv_dz, 0.0f), 1.0f);
    const real old_mass_per_volume = fmaxf(rho_qr_old[local_idx], 0.0f);
    const real old_mass_per_area = old_mass_per_volume * dz;
    const real outgoing_mass_per_area = fallout_fraction * old_mass_per_area;
    const real retained_mass_per_area = old_mass_per_area - outgoing_mass_per_area;
    const real updated_mass_per_area =
        retained_mass_per_area + incoming_mass_per_area;
    rho_qr[local_idx] = updated_mass_per_area * inv_dz;
    incoming_mass_per_area = outgoing_mass_per_area;
  }

  accumulated_surface_precip[static_cast<std::size_t>(j) * nx + i] +=
      incoming_mass_per_area;
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
                                  const WarmRainConfig& config,
                                  std::vector<WarmRainSurfaceAccumulation>*
                                      accumulations) {
  require_warm_rain_tracer_contract(states, tracers);
  if (states.empty()) {
    return;
  }
  if (accumulations != nullptr) {
    gwm::require(accumulations->size() == states.size(),
                 "Warm-rain accumulation/state size mismatch");
    for (std::size_t n = 0; n < states.size(); ++n) {
      if (!(*accumulations)[n].matches(states[n].rho_d.nx(), states[n].rho_d.ny())) {
        (*accumulations)[n].reset(states[n].rho_d.nx(), states[n].rho_d.ny());
      }
    }
  }

  dim3 block(8, 8, 4);
  dim3 column_block(8, 8, 1);
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

    if (accumulations != nullptr && config.rain_terminal_velocity > 0.0f &&
        config.sedimentation_layer_depth_m > 0.0f) {
      dim3 column_grid((state.rho_d.nx() + column_block.x - 1) / column_block.x,
                       (state.rho_d.ny() + column_block.y - 1) / column_block.y,
                       1);
      rain_fallout_nominal_kernel<<<column_grid, column_block>>>(
          tracer.mass(qr_index).data(),
          (*accumulations)[n].liquid_precipitation_kg_m2.data(),
          state.rho_d.nx(), state.rho_d.ny(), state.rho_d.nz(),
          state.rho_d.halo(), config.dt, config.rain_terminal_velocity,
          config.sedimentation_layer_depth_m);
      GWM_CUDA_CHECK(cudaGetLastError());
    }
  }
  GWM_CUDA_CHECK(cudaDeviceSynchronize());
}

void apply_warm_rain_microphysics(
    std::vector<dycore::DryState>& states,
    std::vector<state::TracerState>& tracers,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const WarmRainConfig& config,
    std::vector<WarmRainSurfaceAccumulation>* accumulations) {
  require_warm_rain_tracer_contract(states, tracers);
  gwm::require(states.size() == layout.size(),
               "State/layout size mismatch in warm-rain microphysics");
  gwm::require(!layout.empty(), "Warm-rain microphysics requires nonempty layout");
  gwm::require(metrics.nx == layout.front().global_nx &&
                   metrics.ny == layout.front().global_ny &&
                   metrics.nz == layout.front().nz,
               "GridMetrics/layout mismatch in warm-rain microphysics");
  if (states.empty()) {
    return;
  }
  if (accumulations != nullptr) {
    gwm::require(accumulations->size() == states.size(),
                 "Warm-rain accumulation/state size mismatch");
    for (std::size_t n = 0; n < states.size(); ++n) {
      if (!(*accumulations)[n].matches(states[n].rho_d.nx(), states[n].rho_d.ny())) {
        (*accumulations)[n].reset(states[n].rho_d.nx(), states[n].rho_d.ny(),
                                  metrics.dx * metrics.dy);
      }
    }
  }

  dim3 block(8, 8, 4);
  dim3 column_block(8, 8, 1);
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

    if (accumulations != nullptr && config.rain_terminal_velocity > 0.0f) {
      auto qr_old = tracer.mass(qr_index).clone_empty_like("_qr_old");
      qr_old.copy_all_from(tracer.mass(qr_index));
      dim3 column_grid((state.rho_d.nx() + column_block.x - 1) / column_block.x,
                       (state.rho_d.ny() + column_block.y - 1) / column_block.y,
                       1);
      rain_fallout_metric_kernel<<<column_grid, column_block>>>(
          qr_old.data(), tracer.mass(qr_index).data(),
          (*accumulations)[n].liquid_precipitation_kg_m2.data(),
          state.rho_d.nx(), state.rho_d.ny(), state.rho_d.nz(),
          state.rho_d.halo(), metrics.inv_dz_cell_metric.data(), metrics.nx,
          metrics.ny, layout[n].i_begin, layout[n].j_begin, config.dt,
          config.rain_terminal_velocity);
      GWM_CUDA_CHECK(cudaGetLastError());
    }
  }
  GWM_CUDA_CHECK(cudaDeviceSynchronize());
}

WarmRainSummary summarize_warm_rain(
    const std::vector<state::TracerState>& tracers,
    const std::vector<WarmRainSurfaceAccumulation>* accumulations) {
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
  if (accumulations != nullptr) {
    std::size_t total_cell_count = 0;
    for (const auto& accumulation : *accumulations) {
      for (int j = 0; j < accumulation.ny; ++j) {
        for (int i = 0; i < accumulation.nx; ++i) {
          const double value = static_cast<double>(
              accumulation.liquid_precipitation_kg_m2(i, j, 0));
          summary.accumulated_surface_precipitation_sum_mm += value;
          summary.max_surface_precipitation_mm =
              std::max(summary.max_surface_precipitation_mm, value);
          ++total_cell_count;
        }
      }
    }
    if (total_cell_count > 0) {
      summary.mean_surface_precipitation_mm =
          summary.accumulated_surface_precipitation_sum_mm /
          static_cast<double>(total_cell_count);
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
      << summary.total_rain_water_mass << ",\n";
  oss << indent << "  \"accumulated_surface_precipitation_sum_mm\": "
      << summary.accumulated_surface_precipitation_sum_mm << ",\n";
  oss << indent << "  \"mean_surface_precipitation_mm\": "
      << summary.mean_surface_precipitation_mm << ",\n";
  oss << indent << "  \"max_surface_precipitation_mm\": "
      << summary.max_surface_precipitation_mm << "\n";
  oss << indent << "}";
  return oss.str();
}

}  // namespace gwm::physics
