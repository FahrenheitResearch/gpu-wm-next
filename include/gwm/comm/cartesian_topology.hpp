#pragma once

#include <cstddef>
#include <vector>

#include "gwm/domain/subdomain_descriptor.hpp"

namespace gwm::comm {

struct CartesianNeighborRanks {
  int west = -1;
  int east = -1;
  int south = -1;
  int north = -1;
};

[[nodiscard]] const domain::SubdomainDescriptor* find_rank_descriptor(
    const std::vector<domain::SubdomainDescriptor>& layout, int rank);

[[nodiscard]] CartesianNeighborRanks find_cartesian_neighbors(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::SubdomainDescriptor& desc);

struct ScalarHaloExchangePlan {
  domain::SubdomainDescriptor desc;
  CartesianNeighborRanks neighbors;
  std::size_t x_face_count = 0;
  std::size_t y_face_count = 0;
};

[[nodiscard]] ScalarHaloExchangePlan make_scalar_halo_exchange_plan(
    const domain::SubdomainDescriptor& desc,
    const CartesianNeighborRanks& neighbors);

[[nodiscard]] ScalarHaloExchangePlan make_scalar_halo_exchange_plan(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::SubdomainDescriptor& desc);

}  // namespace gwm::comm
