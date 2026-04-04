#pragma once

#include <vector>

#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/dycore/dry_core.hpp"

namespace gwm::dycore {

[[nodiscard]] bool touches_west_boundary(
    const domain::SubdomainDescriptor& desc);
[[nodiscard]] bool touches_east_boundary(
    const domain::SubdomainDescriptor& desc);
[[nodiscard]] bool touches_south_boundary(
    const domain::SubdomainDescriptor& desc);
[[nodiscard]] bool touches_north_boundary(
    const domain::SubdomainDescriptor& desc);

void apply_reference_boundaries(
    std::vector<DryState>& states, const std::vector<DryState>& boundary_states,
    const std::vector<domain::SubdomainDescriptor>& layout);

}  // namespace gwm::dycore
