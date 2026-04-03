#include "gwm/comm/halo_exchange.hpp"

namespace gwm::comm {

ScalarHaloBuffers HaloExchange::allocate_scalar_buffers(
    const ScalarHaloExchangePlan& plan) {
  ScalarHaloBuffers buffers{};
  buffers.west_send.resize(plan.x_face_count);
  buffers.east_send.resize(plan.x_face_count);
  buffers.south_send.resize(plan.y_face_count);
  buffers.north_send.resize(plan.y_face_count);
  buffers.west_recv.resize(plan.x_face_count);
  buffers.east_recv.resize(plan.x_face_count);
  buffers.south_recv.resize(plan.y_face_count);
  buffers.north_recv.resize(plan.y_face_count);
  return buffers;
}

void HaloExchange::pack_scalar_x_faces(const state::Field3D<real>& field,
                                       const domain::SubdomainDescriptor& desc,
                                       ScalarHaloBuffers& buffers) {
  std::size_t west_idx = 0;
  std::size_t east_idx = 0;
  for (int k = 0; k < desc.nz; ++k) {
    for (int j = 0; j < desc.ny_local(); ++j) {
      for (int h = 0; h < desc.halo; ++h) {
        buffers.west_send[west_idx++] = field(h, j, k);
        buffers.east_send[east_idx++] = field(desc.nx_local() - desc.halo + h, j, k);
      }
    }
  }
}

void HaloExchange::unpack_scalar_x_faces(state::Field3D<real>& field,
                                         const domain::SubdomainDescriptor& desc,
                                         const ScalarHaloBuffers& buffers) {
  std::size_t west_idx = 0;
  std::size_t east_idx = 0;
  for (int k = 0; k < desc.nz; ++k) {
    for (int j = 0; j < desc.ny_local(); ++j) {
      for (int h = 0; h < desc.halo; ++h) {
        field(-desc.halo + h, j, k) = buffers.west_recv[west_idx++];
        field(desc.nx_local() + h, j, k) = buffers.east_recv[east_idx++];
      }
    }
  }
}

void HaloExchange::pack_scalar_y_faces(const state::Field3D<real>& field,
                                       const domain::SubdomainDescriptor& desc,
                                       ScalarHaloBuffers& buffers) {
  std::size_t south_idx = 0;
  std::size_t north_idx = 0;
  for (int k = 0; k < desc.nz; ++k) {
    for (int h = 0; h < desc.halo; ++h) {
      for (int i = -desc.halo; i < desc.nx_local() + desc.halo; ++i) {
        buffers.south_send[south_idx++] = field(i, h, k);
        buffers.north_send[north_idx++] =
            field(i, desc.ny_local() - desc.halo + h, k);
      }
    }
  }
}

void HaloExchange::unpack_scalar_y_faces(state::Field3D<real>& field,
                                         const domain::SubdomainDescriptor& desc,
                                         const ScalarHaloBuffers& buffers) {
  std::size_t south_idx = 0;
  std::size_t north_idx = 0;
  for (int k = 0; k < desc.nz; ++k) {
    for (int h = 0; h < desc.halo; ++h) {
      for (int i = -desc.halo; i < desc.nx_local() + desc.halo; ++i) {
        field(i, -desc.halo + h, k) = buffers.south_recv[south_idx++];
        field(i, desc.ny_local() + h, k) = buffers.north_recv[north_idx++];
      }
    }
  }
}

void HaloExchange::exchange_scalar(
    std::vector<state::Field3D<real>>& fields,
    const std::vector<domain::SubdomainDescriptor>& layout) {
  gwm::require(fields.size() == layout.size(),
               "Field/layout size mismatch in exchange_scalar");

  std::vector<ScalarHaloExchangePlan> plans;
  plans.reserve(layout.size());
  std::vector<ScalarHaloBuffers> buffers;
  buffers.reserve(layout.size());
  for (const auto& desc : layout) {
    plans.push_back(make_scalar_halo_exchange_plan(layout, desc));
    buffers.push_back(allocate_scalar_buffers(plans.back()));
  }

  for (std::size_t n = 0; n < layout.size(); ++n) {
    pack_scalar_x_faces(fields[n], layout[n], buffers[n]);
  }

  for (std::size_t n = 0; n < layout.size(); ++n) {
    const auto& neighbors = plans[n].neighbors;
    if (neighbors.west >= 0) {
      buffers[n].west_recv =
          buffers[static_cast<std::size_t>(neighbors.west)].east_send;
    }
    if (neighbors.east >= 0) {
      buffers[n].east_recv =
          buffers[static_cast<std::size_t>(neighbors.east)].west_send;
    }
    unpack_scalar_x_faces(fields[n], layout[n], buffers[n]);
  }

  for (std::size_t n = 0; n < layout.size(); ++n) {
    pack_scalar_y_faces(fields[n], layout[n], buffers[n]);
  }

  for (std::size_t n = 0; n < layout.size(); ++n) {
    const auto& neighbors = plans[n].neighbors;
    if (neighbors.south >= 0) {
      buffers[n].south_recv =
          buffers[static_cast<std::size_t>(neighbors.south)].north_send;
    }
    if (neighbors.north >= 0) {
      buffers[n].north_recv =
          buffers[static_cast<std::size_t>(neighbors.north)].south_send;
    }
    unpack_scalar_y_faces(fields[n], layout[n], buffers[n]);
  }
}

}  // namespace gwm::comm
