#pragma once

#include <vector>

#include "gwm/comm/virtual_rank_layout.hpp"

namespace gwm::dycore {

struct DryTransportConfig {
  real u_adv = 20.0f;
  real v_adv = 0.0f;
  real dt = 10.0f;
  real dx = 1000.0f;
  real dy = 1000.0f;
};

void advance_scalar_ssprk3(std::vector<state::Field3D<real>>& fields,
                           const std::vector<domain::SubdomainDescriptor>& layout,
                           const DryTransportConfig& config);

}  // namespace gwm::dycore
