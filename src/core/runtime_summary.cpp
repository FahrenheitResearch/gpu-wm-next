#include "gwm/core/runtime_summary.hpp"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>

#include "gwm/core/cuda_utils.hpp"
#include "gwm/state/tracer_registry.hpp"

namespace gwm::core {

namespace {

int tracer_index_or_neg(const std::vector<TracerSummary>& tracers,
                        const std::string& name) {
  for (std::size_t idx = 0; idx < tracers.size(); ++idx) {
    if (tracers[idx].name == name) {
      return static_cast<int>(idx);
    }
  }
  return -1;
}

}  // namespace

RuntimeStateSummary summarize_runtime_state(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics,
    const std::vector<physics::WarmRainSurfaceAccumulation>* accumulations) {
  GWM_CUDA_CHECK(cudaDeviceSynchronize());
  gwm::require(states.size() == tracers.size(),
               "State/tracer size mismatch in summarize_runtime_state");
  gwm::require(states.size() == layout.size(),
               "State/layout size mismatch in summarize_runtime_state");

  RuntimeStateSummary summary{};
  summary.dry = dycore::summarize_dry_states(states);
  if (states.empty() || tracers.empty()) {
    return summary;
  }

  const auto& registry = tracers.front().registry();
  summary.tracers.reserve(registry.size());
  for (std::size_t idx = 0; idx < registry.size(); ++idx) {
    const auto& spec = registry.at(static_cast<int>(idx));
    TracerSummary tracer{};
    tracer.name = spec.name;
    tracer.units = spec.units;
    tracer.positive = spec.positive;
    tracer.column_integrated_sum_kg_m2 = 0.0;
    tracer.mixing_ratio.min = std::numeric_limits<double>::max();
    tracer.mixing_ratio.max = -std::numeric_limits<double>::max();
    summary.tracers.push_back(std::move(tracer));
  }

  bool found_owned_values = false;
  for (std::size_t rank = 0; rank < states.size(); ++rank) {
    gwm::require(tracers[rank].size() == summary.tracers.size(),
                 "Tracer registry mismatch in summarize_runtime_state");
    const auto& desc = layout[rank];
    for (std::size_t idx = 0; idx < summary.tracers.size(); ++idx) {
      const auto& spec = registry.at(static_cast<int>(idx));
      gwm::require(spec.name ==
                       tracers[rank].registry().at(static_cast<int>(idx)).name,
                   "Tracer ordering mismatch in summarize_runtime_state");
      const auto& rho_q = tracers[rank].mass(static_cast<int>(idx));
      for (int k = 0; k < desc.nz; ++k) {
        for (int j = 0; j < desc.ny_local(); ++j) {
          for (int i = 0; i < desc.nx_local(); ++i) {
            const double dz = 1.0 / static_cast<double>(
                metrics.inv_dz_cell(desc.i_begin + i, desc.j_begin + j, k));
            summary.tracers[idx].column_integrated_sum_kg_m2 +=
                static_cast<double>(rho_q(i, j, k)) * dz;
          }
        }
      }
    }

    for (int k = 0; k < states[rank].rho_d.nz(); ++k) {
      for (int j = 0; j < states[rank].rho_d.ny(); ++j) {
        for (int i = 0; i < states[rank].rho_d.nx(); ++i) {
          found_owned_values = true;
          const double rho =
              std::max(static_cast<double>(states[rank].rho_d(i, j, k)), 1.0e-12);
          for (std::size_t idx = 0; idx < summary.tracers.size(); ++idx) {
            const double rho_q =
                static_cast<double>(tracers[rank].mass(static_cast<int>(idx))(i, j, k));
            const double q = rho_q / rho;
            summary.tracers[idx].mixing_ratio.min =
                std::min(summary.tracers[idx].mixing_ratio.min, q);
            summary.tracers[idx].mixing_ratio.max =
                std::max(summary.tracers[idx].mixing_ratio.max, q);
          }
        }
      }
    }
  }

  if (!found_owned_values) {
    for (auto& tracer : summary.tracers) {
      tracer.mixing_ratio.min = 0.0;
      tracer.mixing_ratio.max = 0.0;
    }
  }

  const int qv_index =
      tracer_index_or_neg(summary.tracers, state::kSpecificHumidityTracerName);
  const int qc_index =
      tracer_index_or_neg(summary.tracers, state::kCloudWaterTracerName);
  const int qr_index =
      tracer_index_or_neg(summary.tracers, state::kRainWaterTracerName);

  if (qv_index >= 0) {
    summary.moisture.vapor_water_path_sum_kg_m2 =
        summary.tracers[static_cast<std::size_t>(qv_index)]
            .column_integrated_sum_kg_m2;
  }
  if (qc_index >= 0) {
    summary.moisture.cloud_water_path_sum_kg_m2 =
        summary.tracers[static_cast<std::size_t>(qc_index)]
            .column_integrated_sum_kg_m2;
  }
  if (qr_index >= 0) {
    summary.moisture.rain_water_path_sum_kg_m2 =
        summary.tracers[static_cast<std::size_t>(qr_index)]
            .column_integrated_sum_kg_m2;
  }
  summary.moisture.condensed_water_path_sum_kg_m2 =
      summary.moisture.cloud_water_path_sum_kg_m2 +
      summary.moisture.rain_water_path_sum_kg_m2;
  summary.moisture.total_water_path_sum_kg_m2 =
      summary.moisture.vapor_water_path_sum_kg_m2 +
      summary.moisture.condensed_water_path_sum_kg_m2;
  if (accumulations != nullptr) {
    const auto warm_rain_summary =
        physics::summarize_warm_rain(tracers, layout, metrics, accumulations);
    summary.moisture.accumulated_surface_precipitation_sum_mm =
        warm_rain_summary.accumulated_surface_precipitation_sum_mm;
    summary.moisture.mean_surface_precipitation_mm =
        warm_rain_summary.mean_surface_precipitation_mm;
    summary.moisture.max_surface_precipitation_mm =
        warm_rain_summary.max_surface_precipitation_mm;
    double precipitating_sum_mm = 0.0;
    int total_surface_cell_count = 0;
    int precipitating_surface_cell_count = 0;
    for (const auto& accumulation : *accumulations) {
      for (int j = 0; j < accumulation.ny; ++j) {
        for (int i = 0; i < accumulation.nx; ++i) {
          ++total_surface_cell_count;
          const double value = static_cast<double>(
              accumulation.liquid_precipitation_kg_m2(i, j, 0));
          if (value > 0.0) {
            ++precipitating_surface_cell_count;
            precipitating_sum_mm += value;
          }
        }
      }
    }
    summary.moisture.total_surface_cell_count = total_surface_cell_count;
    summary.moisture.precipitating_surface_cell_count =
        precipitating_surface_cell_count;
    if (total_surface_cell_count > 0) {
      summary.moisture.precipitating_surface_fraction =
          static_cast<double>(precipitating_surface_cell_count) /
          static_cast<double>(total_surface_cell_count);
    }
    if (precipitating_surface_cell_count > 0) {
      summary.moisture.mean_precipitating_surface_precipitation_mm =
          precipitating_sum_mm /
          static_cast<double>(precipitating_surface_cell_count);
    }
  }

  return summary;
}

std::string runtime_state_summary_to_json(const RuntimeStateSummary& summary,
                                          const std::string& indent) {
  std::ostringstream oss;
  oss << std::setprecision(10);
  oss << "{\n";
  oss << indent << "\"total_dry_mass\": " << summary.dry.total_dry_mass << ",\n";
  oss << indent << "\"total_rho_theta_m\": " << summary.dry.total_rho_theta_m
      << ",\n";
  oss << indent << "\"total_horizontal_momentum_x\": "
      << summary.dry.total_horizontal_momentum_x << ",\n";
  oss << indent << "\"total_horizontal_momentum_y\": "
      << summary.dry.total_horizontal_momentum_y << ",\n";
  oss << indent << "\"total_vertical_momentum\": "
      << summary.dry.total_vertical_momentum << ",\n";
  oss << indent << "\"rho_d\": {\"min\": " << summary.dry.rho_d.min
      << ", \"max\": " << summary.dry.rho_d.max << "},\n";
  oss << indent << "\"theta_m\": {\"min\": " << summary.dry.theta_m.min
      << ", \"max\": " << summary.dry.theta_m.max << "},\n";
  oss << indent << "\"u_face\": {\"min\": " << summary.dry.u_face.min
      << ", \"max\": " << summary.dry.u_face.max << "},\n";
  oss << indent << "\"v_face\": {\"min\": " << summary.dry.v_face.min
      << ", \"max\": " << summary.dry.v_face.max << "},\n";
  oss << indent << "\"w_face\": {\"min\": " << summary.dry.w_face.min
      << ", \"max\": " << summary.dry.w_face.max << "},\n";
  oss << indent << "\"moisture\": {\n";
  oss << indent << "  \"vapor_water_path_sum_kg_m2\": "
      << summary.moisture.vapor_water_path_sum_kg_m2 << ",\n";
  oss << indent << "  \"cloud_water_path_sum_kg_m2\": "
      << summary.moisture.cloud_water_path_sum_kg_m2 << ",\n";
  oss << indent << "  \"rain_water_path_sum_kg_m2\": "
      << summary.moisture.rain_water_path_sum_kg_m2 << ",\n";
  oss << indent << "  \"condensed_water_path_sum_kg_m2\": "
      << summary.moisture.condensed_water_path_sum_kg_m2 << ",\n";
  oss << indent << "  \"total_water_path_sum_kg_m2\": "
      << summary.moisture.total_water_path_sum_kg_m2 << ",\n";
  oss << indent << "  \"accumulated_surface_precipitation_sum_mm\": "
      << summary.moisture.accumulated_surface_precipitation_sum_mm << ",\n";
  oss << indent << "  \"mean_surface_precipitation_mm\": "
      << summary.moisture.mean_surface_precipitation_mm << ",\n";
  oss << indent << "  \"max_surface_precipitation_mm\": "
      << summary.moisture.max_surface_precipitation_mm << ",\n";
  oss << indent << "  \"total_surface_cell_count\": "
      << summary.moisture.total_surface_cell_count << ",\n";
  oss << indent << "  \"precipitating_surface_cell_count\": "
      << summary.moisture.precipitating_surface_cell_count << ",\n";
  oss << indent << "  \"precipitating_surface_fraction\": "
      << summary.moisture.precipitating_surface_fraction << ",\n";
  oss << indent << "  \"mean_precipitating_surface_precipitation_mm\": "
      << summary.moisture.mean_precipitating_surface_precipitation_mm << "\n";
  oss << indent << "},\n";
  oss << indent << "\"tracers\": {\n";
  for (std::size_t idx = 0; idx < summary.tracers.size(); ++idx) {
    const auto& tracer = summary.tracers[idx];
    oss << indent << "  \"" << tracer.name << "\": {\n";
    oss << indent << "    \"units\": \"" << tracer.units << "\",\n";
    oss << indent << "    \"positive\": "
        << (tracer.positive ? "true" : "false") << ",\n";
    oss << indent << "    \"column_integrated_sum_kg_m2\": "
        << tracer.column_integrated_sum_kg_m2 << ",\n";
    oss << indent << "    \"mixing_ratio\": {\"min\": "
        << tracer.mixing_ratio.min << ", \"max\": "
        << tracer.mixing_ratio.max << "}\n";
    oss << indent << "  }";
    if (idx + 1 != summary.tracers.size()) {
      oss << ",";
    }
    oss << "\n";
  }
  oss << indent << "}\n";
  oss << "}";
  return oss.str();
}

}  // namespace gwm::core
