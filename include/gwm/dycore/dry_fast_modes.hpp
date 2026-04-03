#pragma once

#include <vector>

#include "gwm/domain/grid_metrics.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/dycore/dry_core.hpp"

namespace gwm::dycore {

void apply_local_split_explicit_fast_modes(
    std::vector<DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const DryStepperConfig& config);

}  // namespace gwm::dycore
