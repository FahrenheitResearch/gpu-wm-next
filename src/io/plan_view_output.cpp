#include "gwm/io/plan_view_output.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/core/cuda_utils.hpp"

namespace gwm::io {

namespace {

using gwm::state::FaceOrientation;

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
          field(i, j, k) = state.rho_theta_m(i, j, k) / std::max(rho, 1.0e-6f);
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
          field(i, j, k) = u_center / std::max(rho, 1.0e-6f);
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
          field(i, j, k) = v_center / std::max(rho, 1.0e-6f);
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
          field(i, j, k) = w_center / std::max(rho, 1.0e-6f);
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
  bundle.fields.push_back(make_slice_field("u_velocity", "m s^-1", u_global, clamped_k));
  bundle.fields.push_back(make_slice_field("v_velocity", "m s^-1", v_global, clamped_k));
  bundle.fields.push_back(
      make_slice_field("wind_speed", "m s^-1", wind_speed_global, clamped_k));
  bundle.fields.push_back(make_slice_field("w_velocity", "m s^-1", w_global, clamped_k));

  return bundle;
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
  std::ofstream out(path, std::ios::binary);
  if (!out) {
    throw std::runtime_error("Failed to open plan-view output path: " + path);
  }
  out << plan_view_bundle_to_json(bundle) << "\n";
}

}  // namespace gwm::io
