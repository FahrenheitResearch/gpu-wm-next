#include "gwm/comm/cartesian_topology.hpp"

#include "gwm/core/types.hpp"

namespace gwm::comm {

namespace {

int wrap_index(int value, int size) {
  int wrapped = value % size;
  if (wrapped < 0) {
    wrapped += size;
  }
  return wrapped;
}

int rank_at(const std::vector<domain::SubdomainDescriptor>& layout, int coord_x,
            int coord_y, bool periodic_x, bool periodic_y) {
  if (layout.empty()) {
    return -1;
  }
  const int ranks_x = layout.front().ranks_x;
  const int ranks_y = layout.front().ranks_y;

  if (periodic_x) {
    coord_x = wrap_index(coord_x, ranks_x);
  } else if (coord_x < 0 || coord_x >= ranks_x) {
    return -1;
  }

  if (periodic_y) {
    coord_y = wrap_index(coord_y, ranks_y);
  } else if (coord_y < 0 || coord_y >= ranks_y) {
    return -1;
  }

  for (const auto& candidate : layout) {
    if (candidate.coord_x == coord_x && candidate.coord_y == coord_y) {
      return candidate.rank;
    }
  }
  return -1;
}

}  // namespace

CartesianNeighborRanks find_cartesian_neighbors(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::SubdomainDescriptor& desc) {
  CartesianNeighborRanks neighbors{};
  neighbors.west =
      rank_at(layout, desc.coord_x - 1, desc.coord_y, desc.periodic_x,
              desc.periodic_y);
  neighbors.east =
      rank_at(layout, desc.coord_x + 1, desc.coord_y, desc.periodic_x,
              desc.periodic_y);
  neighbors.south =
      rank_at(layout, desc.coord_x, desc.coord_y - 1, desc.periodic_x,
              desc.periodic_y);
  neighbors.north =
      rank_at(layout, desc.coord_x, desc.coord_y + 1, desc.periodic_x,
              desc.periodic_y);
  return neighbors;
}

ScalarHaloExchangePlan make_scalar_halo_exchange_plan(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::SubdomainDescriptor& desc) {
  gwm::require(desc.halo >= 0, "Halo width must be nonnegative");

  ScalarHaloExchangePlan plan{};
  plan.desc = desc;
  plan.neighbors = find_cartesian_neighbors(layout, desc);
  plan.x_face_count = static_cast<std::size_t>(desc.halo) * desc.ny_local() *
                      desc.nz;
  plan.y_face_count = static_cast<std::size_t>(desc.halo) *
                      (desc.nx_local() + 2 * desc.halo) * desc.nz;
  return plan;
}

}  // namespace gwm::comm
