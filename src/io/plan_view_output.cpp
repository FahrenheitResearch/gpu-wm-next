#include "gwm/io/plan_view_output.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/core/dry_thermo.hpp"
#include "gwm/core/moist_thermo.hpp"
#include "gwm/core/cuda_utils.hpp"
#include "gwm/ingest/runtime_case.hpp"
#include "gwm/state/tracer_registry.hpp"

namespace gwm::io {

namespace {

using gwm::state::FaceOrientation;

constexpr real kMinRho = 1.0e-6f;
constexpr real kMinPressurePa = 1.0f;

std::vector<state::Field3D<real>> derive_tracer_mixing_ratio_fields(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::string& tracer_name, const std::string& label_suffix);

std::vector<state::Field3D<real>> derive_total_condensate_fields(
    const std::vector<state::Field3D<real>>& cloud_fields,
    const std::vector<state::Field3D<real>>& rain_fields);

std::vector<double> derive_column_rain_water_map(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics);

std::vector<double> derive_column_tracer_map(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const std::string& tracer_name);

std::vector<double> derive_column_total_condensate_map(
    const std::vector<double>& column_cloud_water,
    const std::vector<double>& column_rain_water);

std::vector<double> derive_column_rain_fraction_map(
    const std::vector<double>& column_total_condensate,
    const std::vector<double>& column_rain_water);

std::vector<double> derive_synthetic_reflectivity_map(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics);

std::vector<double> derive_accumulated_surface_precipitation_map(
    const std::vector<physics::WarmRainSurfaceAccumulation>& accumulations,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics);

bool tracer_present_in_all_states(const std::vector<state::TracerState>& tracers,
                                  const std::string& tracer_name) {
  return std::all_of(tracers.begin(), tracers.end(), [&](const auto& tracer_state) {
    return tracer_state.find(tracer_name).has_value();
  });
}

std::vector<state::Field3D<real>> derive_theta_fields(
    const std::vector<dycore::DryState>& states) {
  std::vector<state::Field3D<real>> fields;
  fields.reserve(states.size());
  for (const auto& state : states) {
    auto field = state.rho_d.clone_empty_like("_theta_output");
    for (int k = 0; k < state.rho_d.nz(); ++k) {
      for (int j = 0; j < state.rho_d.ny(); ++j) {
        for (int i = 0; i < state.rho_d.nx(); ++i) {
          const real rho = state.rho_d(i, j, k);
          field(i, j, k) = state.rho_theta_m(i, j, k) / std::max(rho, kMinRho);
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

std::vector<state::Field3D<real>> derive_u_velocity_fields(
    const std::vector<dycore::DryState>& states) {
  std::vector<state::Field3D<real>> fields;
  fields.reserve(states.size());
  for (const auto& state : states) {
    auto field = state.rho_d.clone_empty_like("_u_output");
    const auto& u = state.mom_u.storage();
    for (int k = 0; k < state.rho_d.nz(); ++k) {
      for (int j = 0; j < state.rho_d.ny(); ++j) {
        for (int i = 0; i < state.rho_d.nx(); ++i) {
          const real rho = state.rho_d(i, j, k);
          const real u_center = 0.5f * (u(i, j, k) + u(i + 1, j, k));
          field(i, j, k) = u_center / std::max(rho, kMinRho);
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

std::vector<state::Field3D<real>> derive_v_velocity_fields(
    const std::vector<dycore::DryState>& states) {
  std::vector<state::Field3D<real>> fields;
  fields.reserve(states.size());
  for (const auto& state : states) {
    auto field = state.rho_d.clone_empty_like("_v_output");
    const auto& v = state.mom_v.storage();
    for (int k = 0; k < state.rho_d.nz(); ++k) {
      for (int j = 0; j < state.rho_d.ny(); ++j) {
        for (int i = 0; i < state.rho_d.nx(); ++i) {
          const real rho = state.rho_d(i, j, k);
          const real v_center = 0.5f * (v(i, j, k) + v(i, j + 1, k));
          field(i, j, k) = v_center / std::max(rho, kMinRho);
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

std::vector<state::Field3D<real>> derive_w_velocity_fields(
    const std::vector<dycore::DryState>& states) {
  std::vector<state::Field3D<real>> fields;
  fields.reserve(states.size());
  for (const auto& state : states) {
    auto field = state.rho_d.clone_empty_like("_w_output");
    const auto& w = state.mom_w.storage();
    for (int k = 0; k < state.rho_d.nz(); ++k) {
      for (int j = 0; j < state.rho_d.ny(); ++j) {
        for (int i = 0; i < state.rho_d.nx(); ++i) {
          const real rho = state.rho_d(i, j, k);
          const real w_center = 0.5f * (w(i, j, k) + w(i, j, k + 1));
          field(i, j, k) = w_center / std::max(rho, kMinRho);
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

std::vector<state::Field3D<real>> derive_wind_speed_fields(
    const std::vector<state::Field3D<real>>& u_fields,
    const std::vector<state::Field3D<real>>& v_fields) {
  gwm::require(u_fields.size() == v_fields.size(),
               "u/v field count mismatch in wind-speed derivation");

  std::vector<state::Field3D<real>> fields;
  fields.reserve(u_fields.size());
  for (std::size_t n = 0; n < u_fields.size(); ++n) {
    auto field = u_fields[n].clone_empty_like("_wind_speed_output");
    for (int k = 0; k < field.nz(); ++k) {
      for (int j = 0; j < field.ny(); ++j) {
        for (int i = 0; i < field.nx(); ++i) {
          const real u = u_fields[n](i, j, k);
          const real v = v_fields[n](i, j, k);
          field(i, j, k) = std::sqrt(u * u + v * v);
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

std::vector<state::Field3D<real>> derive_pressure_fields(
    const std::vector<dycore::DryState>& states) {
  std::vector<state::Field3D<real>> fields;
  fields.reserve(states.size());
  for (const auto& state : states) {
    auto field = state.rho_d.clone_empty_like("_pressure_output");
    for (int k = 0; k < state.rho_d.nz(); ++k) {
      for (int j = 0; j < state.rho_d.ny(); ++j) {
        for (int i = 0; i < state.rho_d.nx(); ++i) {
          field(i, j, k) = gwm::core::dry_pressure_from_rho_theta_m(
              state.rho_d(i, j, k), state.rho_theta_m(i, j, k));
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

std::vector<state::Field3D<real>> derive_temperature_fields(
    const std::vector<state::Field3D<real>>& theta_fields,
    const std::vector<state::Field3D<real>>& pressure_fields) {
  gwm::require(theta_fields.size() == pressure_fields.size(),
               "theta/pressure field count mismatch in temperature derivation");

  std::vector<state::Field3D<real>> fields;
  fields.reserve(theta_fields.size());
  for (std::size_t n = 0; n < theta_fields.size(); ++n) {
    auto field = theta_fields[n].clone_empty_like("_temperature_output");
    for (int k = 0; k < field.nz(); ++k) {
      for (int j = 0; j < field.ny(); ++j) {
        for (int i = 0; i < field.nx(); ++i) {
          const real theta = theta_fields[n](i, j, k);
          const real pressure = std::max(pressure_fields[n](i, j, k), kMinPressurePa);
          const real exner = std::pow(
              pressure / gwm::core::kReferencePressure, gwm::core::kKappa);
          field(i, j, k) = theta * exner;
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

std::vector<state::Field3D<real>> derive_specific_humidity_fields(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers) {
  return derive_tracer_mixing_ratio_fields(
      states, tracers, gwm::state::kSpecificHumidityTracerName,
      "_specific_humidity_output");
}

std::vector<state::Field3D<real>> derive_tracer_mixing_ratio_fields(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::string& tracer_name, const std::string& label_suffix) {
  gwm::require(states.size() == tracers.size(),
               "Dry-state/tracer count mismatch in tracer derivation");
  std::vector<state::Field3D<real>> fields;
  fields.reserve(states.size());
  for (std::size_t n = 0; n < states.size(); ++n) {
    const auto tracer_index = tracers[n].find(tracer_name);
    gwm::require(tracer_index.has_value(),
                 tracer_name + " tracer missing from runtime tracer state");
    auto field = states[n].rho_d.clone_empty_like(label_suffix);
    const auto& rho_q = tracers[n].mass(*tracer_index);
    for (int k = 0; k < field.nz(); ++k) {
      for (int j = 0; j < field.ny(); ++j) {
        for (int i = 0; i < field.nx(); ++i) {
          field(i, j, k) = rho_q(i, j, k) /
                           std::max(states[n].rho_d(i, j, k), kMinRho);
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

std::vector<state::Field3D<real>> derive_total_condensate_fields(
    const std::vector<state::Field3D<real>>& cloud_fields,
    const std::vector<state::Field3D<real>>& rain_fields) {
  gwm::require(cloud_fields.size() == rain_fields.size(),
               "cloud/rain field count mismatch in condensate derivation");
  std::vector<state::Field3D<real>> fields;
  fields.reserve(cloud_fields.size());
  for (std::size_t n = 0; n < cloud_fields.size(); ++n) {
    auto field = cloud_fields[n].clone_empty_like("_total_condensate_output");
    for (int k = 0; k < field.nz(); ++k) {
      for (int j = 0; j < field.ny(); ++j) {
        for (int i = 0; i < field.nx(); ++i) {
          field(i, j, k) = cloud_fields[n](i, j, k) + rain_fields[n](i, j, k);
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

std::vector<double> derive_column_rain_water_map(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics) {
  return derive_column_tracer_map(states, tracers, layout, metrics,
                                  gwm::state::kRainWaterTracerName);
}

std::vector<double> derive_column_tracer_map(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const std::string& tracer_name) {
  gwm::require(states.size() == tracers.size() &&
                   states.size() == layout.size(),
               "State/tracer/layout mismatch in column tracer derivation");
  std::vector<std::vector<double>> maps(
      states.size(),
      std::vector<double>(static_cast<std::size_t>(metrics.nx) *
                              static_cast<std::size_t>(metrics.ny),
                          0.0));
  for (std::size_t n = 0; n < states.size(); ++n) {
    const auto tracer_index = tracers[n].find(tracer_name);
    if (!tracer_index.has_value()) {
      continue;
    }
    const auto& rho_q = tracers[n].mass(*tracer_index);
    const auto& desc = layout[n];
    for (int j = 0; j < desc.ny_local(); ++j) {
      for (int i = 0; i < desc.nx_local(); ++i) {
        double column = 0.0;
        for (int k = 0; k < desc.nz; ++k) {
          column += static_cast<double>(rho_q(i, j, k)) *
                    static_cast<double>(
                        1.0f / metrics.inv_dz_cell(desc.i_begin + i,
                                                    desc.j_begin + j, k));
        }
        const auto global_index =
            static_cast<std::size_t>(desc.j_begin + j) *
                static_cast<std::size_t>(metrics.nx) +
            static_cast<std::size_t>(desc.i_begin + i);
        maps[n][global_index] = column;
      }
    }
  }
  std::vector<double> combined(static_cast<std::size_t>(metrics.nx) *
                                   static_cast<std::size_t>(metrics.ny),
                               0.0);
  for (const auto& map : maps) {
    for (std::size_t idx = 0; idx < combined.size(); ++idx) {
      if (map[idx] != 0.0) {
        combined[idx] = map[idx];
      }
    }
  }
  return combined;
}

std::vector<double> derive_column_total_condensate_map(
    const std::vector<double>& column_cloud_water,
    const std::vector<double>& column_rain_water) {
  gwm::require(column_cloud_water.size() == column_rain_water.size(),
               "Column cloud/rain size mismatch in condensate derivation");
  std::vector<double> total(column_cloud_water.size(), 0.0);
  for (std::size_t idx = 0; idx < total.size(); ++idx) {
    total[idx] = column_cloud_water[idx] + column_rain_water[idx];
  }
  return total;
}

std::vector<double> derive_column_rain_fraction_map(
    const std::vector<double>& column_total_condensate,
    const std::vector<double>& column_rain_water) {
  gwm::require(column_total_condensate.size() == column_rain_water.size(),
               "Column condensate/rain size mismatch in rain-fraction derivation");
  std::vector<double> fraction(column_total_condensate.size(), 0.0);
  for (std::size_t idx = 0; idx < fraction.size(); ++idx) {
    const double total = std::max(column_total_condensate[idx], 0.0);
    if (total <= 0.0) {
      continue;
    }
    fraction[idx] =
        std::clamp(column_rain_water[idx] / total, 0.0, 1.0);
  }
  return fraction;
}

std::vector<double> derive_synthetic_reflectivity_map(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics) {
  gwm::require(states.size() == tracers.size() &&
                   states.size() == layout.size(),
               "State/tracer/layout mismatch in reflectivity derivation");
  std::vector<std::vector<double>> maps(
      states.size(),
      std::vector<double>(static_cast<std::size_t>(metrics.nx) *
                              static_cast<std::size_t>(metrics.ny),
                          -20.0));
  for (std::size_t n = 0; n < states.size(); ++n) {
    const auto qr_index = tracers[n].find(gwm::state::kRainWaterTracerName);
    if (!qr_index.has_value()) {
      continue;
    }
    const auto& rho_qr = tracers[n].mass(*qr_index);
    const auto& desc = layout[n];
    for (int j = 0; j < desc.ny_local(); ++j) {
      for (int i = 0; i < desc.nx_local(); ++i) {
        real dbz = -20.0f;
        for (int k = 0; k < desc.nz; ++k) {
          const real rho = states[n].rho_d(i, j, k);
          const real qr = rho_qr(i, j, k) / std::max(rho, kMinRho);
          dbz = std::max(
              dbz, gwm::core::synthetic_reflectivity_dbz_from_rain(rho, qr));
        }
        const auto global_index =
            static_cast<std::size_t>(desc.j_begin + j) *
                static_cast<std::size_t>(metrics.nx) +
            static_cast<std::size_t>(desc.i_begin + i);
        maps[n][global_index] = static_cast<double>(dbz);
      }
    }
  }
  std::vector<double> combined(static_cast<std::size_t>(metrics.nx) *
                                   static_cast<std::size_t>(metrics.ny),
                               -20.0);
  for (const auto& map : maps) {
    for (std::size_t idx = 0; idx < combined.size(); ++idx) {
      combined[idx] = std::max(combined[idx], map[idx]);
    }
  }
  return combined;
}

std::vector<double> derive_accumulated_surface_precipitation_map(
    const std::vector<physics::WarmRainSurfaceAccumulation>& accumulations,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics) {
  gwm::require(accumulations.size() == layout.size(),
               "Accumulation/layout mismatch in surface-precip derivation");
  std::vector<std::vector<double>> maps(
      accumulations.size(),
      std::vector<double>(static_cast<std::size_t>(metrics.nx) *
                              static_cast<std::size_t>(metrics.ny),
                          0.0));
  for (std::size_t n = 0; n < accumulations.size(); ++n) {
    const auto& accumulation = accumulations[n];
    const auto& desc = layout[n];
    gwm::require(accumulation.matches(desc.nx_local(), desc.ny_local()),
                 "Warm-rain accumulation shape mismatch in plan-view output");
    for (int j = 0; j < desc.ny_local(); ++j) {
      for (int i = 0; i < desc.nx_local(); ++i) {
        const auto global_index =
            static_cast<std::size_t>(desc.j_begin + j) *
                static_cast<std::size_t>(metrics.nx) +
            static_cast<std::size_t>(desc.i_begin + i);
        maps[n][global_index] = static_cast<double>(
            accumulation.liquid_precipitation_kg_m2(i, j, 0));
      }
    }
  }
  std::vector<double> combined(static_cast<std::size_t>(metrics.nx) *
                                   static_cast<std::size_t>(metrics.ny),
                               0.0);
  for (const auto& map : maps) {
    for (std::size_t idx = 0; idx < combined.size(); ++idx) {
      if (map[idx] != 0.0) {
        combined[idx] = map[idx];
      }
    }
  }
  return combined;
}

std::vector<state::Field3D<real>> derive_relative_humidity_fields(
    const std::vector<state::Field3D<real>>& q_fields,
    const std::vector<state::Field3D<real>>& pressure_fields,
    const std::vector<state::Field3D<real>>& temperature_fields) {
  gwm::require(q_fields.size() == pressure_fields.size() &&
                   q_fields.size() == temperature_fields.size(),
               "q/pressure/temperature field count mismatch in RH derivation");
  std::vector<state::Field3D<real>> fields;
  fields.reserve(q_fields.size());
  for (std::size_t n = 0; n < q_fields.size(); ++n) {
    auto field = q_fields[n].clone_empty_like("_relative_humidity_output");
    for (int k = 0; k < field.nz(); ++k) {
      for (int j = 0; j < field.ny(); ++j) {
        for (int i = 0; i < field.nx(); ++i) {
          field(i, j, k) = core::relative_humidity_from_specific_humidity(
              q_fields[n](i, j, k), temperature_fields[n](i, j, k),
              pressure_fields[n](i, j, k));
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

std::vector<state::Field3D<real>> derive_dewpoint_fields(
    const std::vector<state::Field3D<real>>& q_fields,
    const std::vector<state::Field3D<real>>& pressure_fields) {
  gwm::require(q_fields.size() == pressure_fields.size(),
               "q/pressure field count mismatch in dewpoint derivation");
  std::vector<state::Field3D<real>> fields;
  fields.reserve(q_fields.size());
  for (std::size_t n = 0; n < q_fields.size(); ++n) {
    auto field = q_fields[n].clone_empty_like("_dewpoint_output");
    for (int k = 0; k < field.nz(); ++k) {
      for (int j = 0; j < field.ny(); ++j) {
        for (int i = 0; i < field.nx(); ++i) {
          field(i, j, k) = core::dewpoint_from_specific_humidity(
              q_fields[n](i, j, k), pressure_fields[n](i, j, k));
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

PlanViewField make_slice_field(const std::string& name, const std::string& units,
                               const state::Field3D<real>& global_field,
                               int slice_k) {
  PlanViewField out{};
  out.name = name;
  out.units = units;
  out.location = "cell_center";
  out.nx = global_field.nx();
  out.ny = global_field.ny();
  out.values.reserve(static_cast<std::size_t>(out.nx) *
                     static_cast<std::size_t>(out.ny));
  for (int j = 0; j < out.ny; ++j) {
    for (int i = 0; i < out.nx; ++i) {
      out.values.push_back(static_cast<double>(global_field(i, j, slice_k)));
    }
  }
  return out;
}

PlanViewField make_flat_slice_field(const std::string& name,
                                    const std::string& units, int nx, int ny,
                                    std::vector<double> values) {
  gwm::require(static_cast<int>(values.size()) == nx * ny,
               "Plan-view slice field size mismatch for " + name);
  PlanViewField out{};
  out.name = name;
  out.units = units;
  out.location = "cell_center";
  out.nx = nx;
  out.ny = ny;
  out.values = std::move(values);
  return out;
}

PlanViewField make_metric_slice_field(const std::string& name,
                                      const std::string& units,
                                      const domain::GridMetrics& metrics,
                                      int slice_k) {
  PlanViewField out{};
  out.name = name;
  out.units = units;
  out.location = "cell_center";
  out.nx = metrics.nx;
  out.ny = metrics.ny;
  out.values.reserve(static_cast<std::size_t>(out.nx) *
                     static_cast<std::size_t>(out.ny));
  for (int j = 0; j < out.ny; ++j) {
    for (int i = 0; i < out.nx; ++i) {
      out.values.push_back(static_cast<double>(metrics.z_center(i, j, slice_k)));
    }
  }
  return out;
}

PlanViewField make_terrain_field(const domain::GridMetrics& metrics) {
  PlanViewField out{};
  out.name = "terrain_height";
  out.units = "m";
  out.location = "cell_center";
  out.nx = metrics.nx;
  out.ny = metrics.ny;
  out.values.reserve(static_cast<std::size_t>(out.nx) *
                     static_cast<std::size_t>(out.ny));
  for (int j = 0; j < out.ny; ++j) {
    for (int i = 0; i < out.nx; ++i) {
      out.values.push_back(static_cast<double>(metrics.terrain_height(i, j)));
    }
  }
  return out;
}

const PlanViewField* find_field(const PlanViewBundle& bundle,
                                const std::string& name) {
  const auto it = std::find_if(bundle.fields.begin(), bundle.fields.end(),
                               [&](const auto& field) {
                                 return field.name == name;
                               });
  return it == bundle.fields.end() ? nullptr : &(*it);
}

void set_or_append_field(PlanViewBundle& bundle, PlanViewField field) {
  const auto it = std::find_if(bundle.fields.begin(), bundle.fields.end(),
                               [&](const auto& existing) {
                                 return existing.name == field.name;
                               });
  if (it != bundle.fields.end()) {
    *it = std::move(field);
  } else {
    bundle.fields.push_back(std::move(field));
  }
}

std::size_t linear_index_3d(int i, int j, int k, int nx, int ny) {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

void maybe_enrich_bundle_from_companion_analysis(PlanViewBundle& bundle,
                                                 const std::string& output_path) {
  const auto parent = std::filesystem::path(output_path).parent_path();
  if (parent.empty()) {
    return;
  }
  const auto analysis_path = parent / "analysis_state.json";
  if (!std::filesystem::exists(analysis_path)) {
    return;
  }
  enrich_plan_view_bundle_from_prepared_case(bundle, analysis_path.string());
}

void append_json_string(std::ostringstream& oss, const std::string& value) {
  oss << "\"";
  for (const char ch : value) {
    if (ch == '\\' || ch == '"') {
      oss << '\\';
    }
    oss << ch;
  }
  oss << "\"";
}

}  // namespace

PlanViewBundle extract_dry_plan_view(
    const std::vector<dycore::DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const std::string& case_kind, int steps,
    real dt, int slice_k) {
  gwm::require(!states.empty(), "extract_dry_plan_view requires at least one state");
  gwm::require(states.size() == layout.size(),
               "State/layout size mismatch in extract_dry_plan_view");

  GWM_CUDA_CHECK(cudaDeviceSynchronize());

  const int clamped_k = std::clamp(slice_k, 0, metrics.nz - 1);

  const auto theta_fields = derive_theta_fields(states);
  const auto pressure_fields = derive_pressure_fields(states);
  const auto temperature_fields =
      derive_temperature_fields(theta_fields, pressure_fields);
  const auto u_fields = derive_u_velocity_fields(states);
  const auto v_fields = derive_v_velocity_fields(states);
  const auto w_fields = derive_w_velocity_fields(states);
  const auto wind_speed_fields = derive_wind_speed_fields(u_fields, v_fields);

  const auto rho_global =
      comm::VirtualRankLayout::gather_scalar([&]() {
        std::vector<state::Field3D<real>> locals;
        locals.reserve(states.size());
        for (const auto& state : states) {
          auto field = state.rho_d.clone_empty_like("_rho_output_copy");
          field.copy_all_from(state.rho_d);
          locals.push_back(std::move(field));
        }
        return locals;
      }(),
      layout, "rho_plan_view");
  const auto theta_global =
      comm::VirtualRankLayout::gather_scalar(theta_fields, layout, "theta_plan_view");
  const auto pressure_global = comm::VirtualRankLayout::gather_scalar(
      pressure_fields, layout, "pressure_plan_view");
  const auto temperature_global = comm::VirtualRankLayout::gather_scalar(
      temperature_fields, layout, "temperature_plan_view");
  const auto u_global =
      comm::VirtualRankLayout::gather_scalar(u_fields, layout, "u_plan_view");
  const auto v_global =
      comm::VirtualRankLayout::gather_scalar(v_fields, layout, "v_plan_view");
  const auto w_global =
      comm::VirtualRankLayout::gather_scalar(w_fields, layout, "w_plan_view");
  const auto wind_speed_global = comm::VirtualRankLayout::gather_scalar(
      wind_speed_fields, layout, "wind_speed_plan_view");

  PlanViewBundle bundle{};
  bundle.case_kind = case_kind;
  bundle.steps = steps;
  bundle.dt = dt;
  bundle.slice_k = clamped_k;
  bundle.nx = metrics.nx;
  bundle.ny = metrics.ny;
  bundle.dx = metrics.dx;
  bundle.dy = metrics.dy;

  double slice_height_sum = 0.0;
  for (int j = 0; j < metrics.ny; ++j) {
    for (int i = 0; i < metrics.nx; ++i) {
      slice_height_sum += static_cast<double>(metrics.z_center(i, j, clamped_k));
    }
  }
  bundle.slice_mean_height_m =
      slice_height_sum / static_cast<double>(metrics.nx * metrics.ny);

  bundle.fields.push_back(make_terrain_field(metrics));
  bundle.fields.push_back(make_metric_slice_field("z_center", "m", metrics, clamped_k));
  bundle.fields.push_back(make_slice_field("rho_d", "kg m^-3", rho_global, clamped_k));
  bundle.fields.push_back(make_slice_field("theta_m", "K", theta_global, clamped_k));
  bundle.fields.push_back(
      make_slice_field("air_pressure", "Pa", pressure_global, clamped_k));
  bundle.fields.push_back(
      make_slice_field("air_temperature", "K", temperature_global, clamped_k));
  bundle.fields.push_back(make_slice_field("u_velocity", "m s^-1", u_global, clamped_k));
  bundle.fields.push_back(make_slice_field("v_velocity", "m s^-1", v_global, clamped_k));
  bundle.fields.push_back(
      make_slice_field("wind_speed", "m s^-1", wind_speed_global, clamped_k));
  bundle.fields.push_back(make_slice_field("w_velocity", "m s^-1", w_global, clamped_k));

  return bundle;
}

PlanViewBundle extract_runtime_plan_view(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const std::string& case_kind, int steps,
    real dt, int slice_k,
    const std::vector<physics::WarmRainSurfaceAccumulation>* accumulations) {
  auto bundle =
      extract_dry_plan_view(states, layout, metrics, case_kind, steps, dt, slice_k);
  if (tracers.empty()) {
    return bundle;
  }

  const int clamped_k = std::clamp(slice_k, 0, metrics.nz - 1);
  const auto pressure_fields = derive_pressure_fields(states);
  const auto temperature_fields =
      derive_temperature_fields(derive_theta_fields(states), pressure_fields);
  const auto q_fields = derive_specific_humidity_fields(states, tracers);
  const auto rh_fields =
      derive_relative_humidity_fields(q_fields, pressure_fields, temperature_fields);
  const auto dewpoint_fields = derive_dewpoint_fields(q_fields, pressure_fields);

  const auto q_global =
      comm::VirtualRankLayout::gather_scalar(q_fields, layout, "specific_humidity_plan_view");
  const auto rh_global = comm::VirtualRankLayout::gather_scalar(
      rh_fields, layout, "relative_humidity_plan_view");
  const auto dewpoint_global = comm::VirtualRankLayout::gather_scalar(
      dewpoint_fields, layout, "dewpoint_plan_view");

  set_or_append_field(bundle, make_slice_field("specific_humidity", "kg kg^-1",
                                               q_global, clamped_k));
  set_or_append_field(bundle,
                      make_slice_field("relative_humidity", "%", rh_global,
                                       clamped_k));
  set_or_append_field(bundle,
                      make_slice_field("dewpoint", "K", dewpoint_global,
                                       clamped_k));
  if (tracer_present_in_all_states(tracers, gwm::state::kCloudWaterTracerName) &&
      tracer_present_in_all_states(tracers, gwm::state::kRainWaterTracerName)) {
    const auto qc_fields = derive_tracer_mixing_ratio_fields(
        states, tracers, gwm::state::kCloudWaterTracerName,
        "_cloud_water_output");
    const auto qr_fields = derive_tracer_mixing_ratio_fields(
        states, tracers, gwm::state::kRainWaterTracerName,
        "_rain_water_output");
    const auto condensate_fields =
        derive_total_condensate_fields(qc_fields, qr_fields);
    const auto qc_global = comm::VirtualRankLayout::gather_scalar(
        qc_fields, layout, "cloud_water_plan_view");
    const auto qr_global = comm::VirtualRankLayout::gather_scalar(
        qr_fields, layout, "rain_water_plan_view");
    const auto condensate_global = comm::VirtualRankLayout::gather_scalar(
        condensate_fields, layout, "condensate_plan_view");
    const auto column_cloud = derive_column_tracer_map(
        states, tracers, layout, metrics, gwm::state::kCloudWaterTracerName);
    const auto column_rain = derive_column_rain_water_map(
        states, tracers, layout, metrics);
    const auto column_total_condensate =
        derive_column_total_condensate_map(column_cloud, column_rain);
    const auto column_rain_fraction = derive_column_rain_fraction_map(
        column_total_condensate, column_rain);
    const auto reflectivity = derive_synthetic_reflectivity_map(
        states, tracers, layout, metrics);

    set_or_append_field(bundle,
                        make_slice_field(gwm::state::kCloudWaterTracerName,
                                         "kg kg^-1", qc_global, clamped_k));
    set_or_append_field(bundle,
                        make_slice_field(gwm::state::kRainWaterTracerName,
                                         "kg kg^-1", qr_global, clamped_k));
    set_or_append_field(bundle,
                        make_slice_field("total_condensate", "kg kg^-1",
                                         condensate_global, clamped_k));
    set_or_append_field(bundle,
                        make_flat_slice_field("column_cloud_water", "kg m^-2",
                                              metrics.nx, metrics.ny,
                                              std::move(column_cloud)));
    set_or_append_field(bundle,
                        make_flat_slice_field("column_rain_water", "kg m^-2",
                                              metrics.nx, metrics.ny,
                                              std::move(column_rain)));
    set_or_append_field(bundle,
                        make_flat_slice_field("column_total_condensate",
                                              "kg m^-2", metrics.nx,
                                              metrics.ny,
                                              std::move(column_total_condensate)));
    set_or_append_field(bundle,
                        make_flat_slice_field("column_rain_fraction", "1",
                                              metrics.nx, metrics.ny,
                                              std::move(column_rain_fraction)));
    set_or_append_field(bundle,
                        make_flat_slice_field("synthetic_reflectivity", "dBZ",
                                              metrics.nx, metrics.ny,
                                              std::move(reflectivity)));
    if (accumulations != nullptr) {
      const auto accumulated_surface_precip =
          derive_accumulated_surface_precipitation_map(*accumulations, layout,
                                                       metrics);
      std::vector<double> mean_surface_precipitation_rate;
      mean_surface_precipitation_rate.reserve(accumulated_surface_precip.size());
      const double elapsed_hours =
          std::max(static_cast<double>(steps) * static_cast<double>(dt) / 3600.0,
                   1.0e-9);
      for (const double value : accumulated_surface_precip) {
        mean_surface_precipitation_rate.push_back(value / elapsed_hours);
      }
      set_or_append_field(bundle,
                          make_flat_slice_field(
                              "accumulated_surface_precipitation", "mm",
                              metrics.nx, metrics.ny,
                              std::move(accumulated_surface_precip)));
      set_or_append_field(bundle,
                          make_flat_slice_field(
                              "mean_surface_precipitation_rate", "mm h^-1",
                              metrics.nx, metrics.ny,
                              std::move(mean_surface_precipitation_rate)));
    }
  }
  return bundle;
}

void enrich_plan_view_bundle_from_prepared_case(
    PlanViewBundle& bundle, const std::string& analysis_state_path) {
  const auto analysis = ingest::load_analysis_state_json(analysis_state_path);
  gwm::require(analysis.grid.nx == bundle.nx && analysis.grid.ny == bundle.ny,
               "Analysis grid mismatch in plan-view enrichment");
  gwm::require(bundle.slice_k >= 0 && bundle.slice_k < analysis.grid.nz,
               "Plan-view slice_k is outside analysis grid for enrichment");

  const auto* pressure_field = find_field(bundle, "air_pressure");
  const auto* temperature_field = find_field(bundle, "air_temperature");
  gwm::require(pressure_field != nullptr && temperature_field != nullptr,
               "Plan-view enrichment requires air_pressure and air_temperature");

  if (find_field(bundle, "specific_humidity") == nullptr) {
    const auto q_it = analysis.atmosphere.values.find("specific_humidity");
    gwm::require(q_it != analysis.atmosphere.values.end(),
                 "Prepared-case analysis is missing specific_humidity");

    std::vector<double> q_values;
    std::vector<double> rh_values;
    std::vector<double> dewpoint_values;
    const auto cell_count = static_cast<std::size_t>(analysis.grid.nx) *
                            static_cast<std::size_t>(analysis.grid.ny);
    q_values.reserve(cell_count);
    rh_values.reserve(cell_count);
    dewpoint_values.reserve(cell_count);

    for (int j = 0; j < analysis.grid.ny; ++j) {
      for (int i = 0; i < analysis.grid.nx; ++i) {
        const auto idx = linear_index_3d(i, j, bundle.slice_k, analysis.grid.nx,
                                         analysis.grid.ny);
        const real q = std::max(q_it->second[idx], 0.0f);
        const auto flat_idx =
            static_cast<std::size_t>(j) * static_cast<std::size_t>(analysis.grid.nx) +
            static_cast<std::size_t>(i);
        const real pressure = static_cast<real>(pressure_field->values.at(flat_idx));
        const real temperature =
            static_cast<real>(temperature_field->values.at(flat_idx));
        q_values.push_back(static_cast<double>(q));
        rh_values.push_back(static_cast<double>(
            core::relative_humidity_from_specific_humidity(
                q, temperature, std::max(pressure, kMinPressurePa))));
        dewpoint_values.push_back(static_cast<double>(
            core::dewpoint_from_specific_humidity(
                q, std::max(pressure, kMinPressurePa))));
      }
    }

    set_or_append_field(bundle,
                        make_flat_slice_field("specific_humidity", "kg kg^-1",
                                              bundle.nx, bundle.ny,
                                              std::move(q_values)));
    set_or_append_field(bundle,
                        make_flat_slice_field("relative_humidity", "%",
                                              bundle.nx, bundle.ny,
                                              std::move(rh_values)));
    set_or_append_field(bundle,
                        make_flat_slice_field("dewpoint", "K", bundle.nx,
                                              bundle.ny,
                                              std::move(dewpoint_values)));
  }

  const auto& surface_fields = analysis.surface.values;
  if (const auto it = surface_fields.find("surface_pressure");
      it != surface_fields.end()) {
    set_or_append_field(bundle,
                        make_flat_slice_field("surface_pressure", "Pa",
                                              bundle.nx, bundle.ny,
                                              std::vector<double>(it->second.begin(),
                                                                  it->second.end())));
  }
  if (const auto it = surface_fields.find("air_temperature_2m");
      it != surface_fields.end()) {
    set_or_append_field(bundle,
                        make_flat_slice_field("air_temperature_2m", "K",
                                              bundle.nx, bundle.ny,
                                              std::vector<double>(it->second.begin(),
                                                                  it->second.end())));
  }
  if (const auto q2_it = surface_fields.find("specific_humidity_2m");
      q2_it != surface_fields.end()) {
    set_or_append_field(bundle,
                        make_flat_slice_field("specific_humidity_2m",
                                              "kg kg^-1", bundle.nx, bundle.ny,
                                              std::vector<double>(q2_it->second.begin(),
                                                                  q2_it->second.end())));
    if (const auto p2_it = surface_fields.find("surface_pressure");
        p2_it != surface_fields.end()) {
      const auto& t2 = surface_fields.at("air_temperature_2m");
      std::vector<double> rh2_values;
      std::vector<double> dewpoint2_values;
      rh2_values.reserve(q2_it->second.size());
      dewpoint2_values.reserve(q2_it->second.size());
      for (std::size_t idx = 0; idx < q2_it->second.size(); ++idx) {
        const real q2 = std::max(q2_it->second[idx], 0.0f);
        const real pressure = std::max(p2_it->second[idx], kMinPressurePa);
        rh2_values.push_back(static_cast<double>(
            core::relative_humidity_from_specific_humidity(
                q2, t2[idx], pressure)));
        dewpoint2_values.push_back(static_cast<double>(
            core::dewpoint_from_specific_humidity(q2, pressure)));
      }
      set_or_append_field(bundle,
                          make_flat_slice_field("relative_humidity_2m", "%",
                                                bundle.nx, bundle.ny,
                                                std::move(rh2_values)));
      set_or_append_field(bundle,
                          make_flat_slice_field("dewpoint_2m", "K", bundle.nx,
                                                bundle.ny,
                                                std::move(dewpoint2_values)));
    }
  }
}

std::string plan_view_bundle_to_json(const PlanViewBundle& bundle) {
  std::ostringstream oss;
  oss << std::setprecision(10);
  oss << "{\n";
  oss << "  \"schema_version\": ";
  append_json_string(oss, bundle.schema_version);
  oss << ",\n";
  oss << "  \"case\": ";
  append_json_string(oss, bundle.case_kind);
  oss << ",\n";
  oss << "  \"steps\": " << bundle.steps << ",\n";
  oss << "  \"dt\": " << bundle.dt << ",\n";
  oss << "  \"slice_k\": " << bundle.slice_k << ",\n";
  oss << "  \"grid\": {\n";
  oss << "    \"nx\": " << bundle.nx << ",\n";
  oss << "    \"ny\": " << bundle.ny << ",\n";
  oss << "    \"dx\": " << bundle.dx << ",\n";
  oss << "    \"dy\": " << bundle.dy << ",\n";
  oss << "    \"slice_mean_height_m\": " << bundle.slice_mean_height_m << "\n";
  oss << "  },\n";
  oss << "  \"fields\": [\n";
  for (std::size_t n = 0; n < bundle.fields.size(); ++n) {
    const auto& field = bundle.fields[n];
    oss << "    {\n";
    oss << "      \"name\": ";
    append_json_string(oss, field.name);
    oss << ",\n";
    oss << "      \"units\": ";
    append_json_string(oss, field.units);
    oss << ",\n";
    oss << "      \"location\": ";
    append_json_string(oss, field.location);
    oss << ",\n";
    oss << "      \"nx\": " << field.nx << ",\n";
    oss << "      \"ny\": " << field.ny << ",\n";
    oss << "      \"storage\": \"row_major_yx\",\n";
    oss << "      \"values\": [";
    for (std::size_t idx = 0; idx < field.values.size(); ++idx) {
      if (idx > 0) {
        oss << ", ";
      }
      oss << field.values[idx];
    }
    oss << "]\n";
    oss << "    }";
    if (n + 1 < bundle.fields.size()) {
      oss << ",";
    }
    oss << "\n";
  }
  oss << "  ]\n";
  oss << "}";
  return oss.str();
}

void write_plan_view_bundle_json(const PlanViewBundle& bundle,
                                 const std::string& path) {
  PlanViewBundle output_bundle = bundle;
  maybe_enrich_bundle_from_companion_analysis(output_bundle, path);
  std::ofstream out(path, std::ios::binary);
  if (!out) {
    throw std::runtime_error("Failed to open plan-view output path: " + path);
  }
  out << plan_view_bundle_to_json(output_bundle) << "\n";
}

}  // namespace gwm::io
