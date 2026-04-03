#pragma once

#include <functional>
#include <vector>

#include "gwm/core/types.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/state/field3d.hpp"

namespace gwm::comm {

class VirtualRankLayout {
 public:
  static std::vector<domain::SubdomainDescriptor> build(
      int global_nx, int global_ny, int nz, int halo, int ranks_x, int ranks_y,
      bool periodic_x = true, bool periodic_y = true);

  static std::vector<state::Field3D<real>> scatter_scalar(
      const std::vector<domain::SubdomainDescriptor>& layout,
      const std::function<real(int, int, int)>& initializer,
      const std::string& label_prefix);

  static state::Field3D<real> gather_scalar(
      const std::vector<state::Field3D<real>>& fields,
      const std::vector<domain::SubdomainDescriptor>& layout,
      const std::string& label = "gathered");
};

}  // namespace gwm::comm
