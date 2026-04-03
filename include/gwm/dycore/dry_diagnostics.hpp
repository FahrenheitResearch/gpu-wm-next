#pragma once

#include <string>
#include <vector>

#include "gwm/dycore/dry_core.hpp"

namespace gwm::dycore {

struct ScalarRange {
  double min = 0.0;
  double max = 0.0;
};

struct DryStateSummary {
  double total_dry_mass = 0.0;
  double total_rho_theta_m = 0.0;
  double total_vertical_momentum = 0.0;
  ScalarRange rho_d;
  ScalarRange theta_m;
  ScalarRange w_face;
};

[[nodiscard]] DryStateSummary summarize_dry_states(
    const std::vector<DryState>& states);

[[nodiscard]] std::string summary_to_json(const DryStateSummary& summary,
                                          const std::string& indent = "  ");

}  // namespace gwm::dycore
