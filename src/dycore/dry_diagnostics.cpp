#include "gwm/dycore/dry_diagnostics.hpp"

#include <algorithm>
#include <limits>
#include <sstream>

#include "gwm/core/cuda_utils.hpp"

namespace gwm::dycore {

DryStateSummary summarize_dry_states(const std::vector<DryState>& states) {
  GWM_CUDA_CHECK(cudaDeviceSynchronize());

  DryStateSummary summary{};
  summary.rho_d.min = std::numeric_limits<double>::max();
  summary.rho_d.max = -std::numeric_limits<double>::max();
  summary.theta_m.min = std::numeric_limits<double>::max();
  summary.theta_m.max = -std::numeric_limits<double>::max();
  summary.u_face.min = std::numeric_limits<double>::max();
  summary.u_face.max = -std::numeric_limits<double>::max();
  summary.v_face.min = std::numeric_limits<double>::max();
  summary.v_face.max = -std::numeric_limits<double>::max();
  summary.w_face.min = std::numeric_limits<double>::max();
  summary.w_face.max = -std::numeric_limits<double>::max();

  for (const auto& state : states) {
    summary.total_dry_mass += state.total_dry_mass();
    summary.total_rho_theta_m += state.total_rho_theta_m();
    summary.total_horizontal_momentum_x += state.mom_u.storage().owned_sum();
    summary.total_horizontal_momentum_y += state.mom_v.storage().owned_sum();
    summary.total_vertical_momentum += state.total_vertical_momentum();

    for (int k = 0; k < state.rho_d.nz(); ++k) {
      for (int j = 0; j < state.rho_d.ny(); ++j) {
        for (int i = 0; i < state.rho_d.nx(); ++i) {
          const double rho = static_cast<double>(state.rho_d(i, j, k));
          const double theta =
              static_cast<double>(state.rho_theta_m(i, j, k)) /
              std::max(rho, 1.0e-12);
          summary.rho_d.min = std::min(summary.rho_d.min, rho);
          summary.rho_d.max = std::max(summary.rho_d.max, rho);
          summary.theta_m.min = std::min(summary.theta_m.min, theta);
          summary.theta_m.max = std::max(summary.theta_m.max, theta);
        }
      }
    }

    const auto& u = state.mom_u.storage();
    for (int k = 0; k < u.nz(); ++k) {
      for (int j = 0; j < u.ny(); ++j) {
        for (int i = 0; i < u.nx(); ++i) {
          const double value = static_cast<double>(u(i, j, k));
          summary.u_face.min = std::min(summary.u_face.min, value);
          summary.u_face.max = std::max(summary.u_face.max, value);
        }
      }
    }

    const auto& v = state.mom_v.storage();
    for (int k = 0; k < v.nz(); ++k) {
      for (int j = 0; j < v.ny(); ++j) {
        for (int i = 0; i < v.nx(); ++i) {
          const double value = static_cast<double>(v(i, j, k));
          summary.v_face.min = std::min(summary.v_face.min, value);
          summary.v_face.max = std::max(summary.v_face.max, value);
        }
      }
    }

    const auto& w = state.mom_w.storage();
    for (int k = 0; k < w.nz(); ++k) {
      for (int j = 0; j < w.ny(); ++j) {
        for (int i = 0; i < w.nx(); ++i) {
          const double value = static_cast<double>(w(i, j, k));
          summary.w_face.min = std::min(summary.w_face.min, value);
          summary.w_face.max = std::max(summary.w_face.max, value);
        }
      }
    }
  }

  if (states.empty()) {
    summary.rho_d.min = 0.0;
    summary.rho_d.max = 0.0;
    summary.theta_m.min = 0.0;
    summary.theta_m.max = 0.0;
    summary.u_face.min = 0.0;
    summary.u_face.max = 0.0;
    summary.v_face.min = 0.0;
    summary.v_face.max = 0.0;
    summary.w_face.min = 0.0;
    summary.w_face.max = 0.0;
  }

  return summary;
}

std::string summary_to_json(const DryStateSummary& summary,
                            const std::string& indent) {
  std::ostringstream oss;
  oss << "{\n";
  oss << indent << "\"total_dry_mass\": " << summary.total_dry_mass << ",\n";
  oss << indent << "\"total_rho_theta_m\": " << summary.total_rho_theta_m
      << ",\n";
  oss << indent << "\"total_horizontal_momentum_x\": "
      << summary.total_horizontal_momentum_x << ",\n";
  oss << indent << "\"total_horizontal_momentum_y\": "
      << summary.total_horizontal_momentum_y << ",\n";
  oss << indent << "\"total_vertical_momentum\": "
      << summary.total_vertical_momentum << ",\n";
  oss << indent << "\"rho_d\": {\"min\": " << summary.rho_d.min
      << ", \"max\": " << summary.rho_d.max << "},\n";
  oss << indent << "\"theta_m\": {\"min\": " << summary.theta_m.min
      << ", \"max\": " << summary.theta_m.max << "},\n";
  oss << indent << "\"u_face\": {\"min\": " << summary.u_face.min
      << ", \"max\": " << summary.u_face.max << "},\n";
  oss << indent << "\"v_face\": {\"min\": " << summary.v_face.min
      << ", \"max\": " << summary.v_face.max << "},\n";
  oss << indent << "\"w_face\": {\"min\": " << summary.w_face.min
      << ", \"max\": " << summary.w_face.max << "}\n";
  oss << "}";
  return oss.str();
}

}  // namespace gwm::dycore
