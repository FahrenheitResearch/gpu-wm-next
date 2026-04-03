#pragma once

#include <vector>

#include "gwm/domain/grid_metrics.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/state/field3d.hpp"

namespace gwm::dycore {

void compute_dry_pressure_fields(
    const std::vector<state::Field3D<real>>& rho_fields,
    const std::vector<state::Field3D<real>>& rho_theta_fields,
    std::vector<state::Field3D<real>>& pressure_fields);

void add_horizontal_pressure_gradient_tendencies(
    const std::vector<state::Field3D<real>>& pressure_fields,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, std::vector<DrySlowTendencies>& out);

}  // namespace gwm::dycore
