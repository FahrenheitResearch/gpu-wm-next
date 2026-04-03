#pragma once

#include <vector>

#include "gwm/domain/grid_metrics.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/state/face_field.hpp"
#include "gwm/state/field3d.hpp"

namespace gwm::dycore {

void add_dry_momentum_flux_tendencies(
    const std::vector<state::Field3D<real>>& rho_fields,
    const std::vector<state::FaceField<real>>& mom_u_fields,
    const std::vector<state::FaceField<real>>& mom_v_fields,
    const std::vector<state::FaceField<real>>& mom_w_fields,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, std::vector<DrySlowTendencies>& out);

}  // namespace gwm::dycore
