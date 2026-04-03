#include "gwm/comm/halo_exchange.hpp"

namespace gwm::comm {

namespace {

struct FaceHaloShape {
  int owned_nx = 0;
  int owned_ny = 0;
  int owned_nz = 0;
  std::size_t x_face_count = 0;
  std::size_t y_face_count = 0;
};

FaceHaloShape make_face_halo_shape(const state::FaceField<real>& field,
                                   const domain::SubdomainDescriptor& desc) {
  FaceHaloShape shape{};
  shape.owned_nx = field.storage().nx();
  shape.owned_ny = field.storage().ny();
  shape.owned_nz = field.storage().nz();
  shape.x_face_count = static_cast<std::size_t>(desc.halo) * shape.owned_ny *
                       shape.owned_nz;
  shape.y_face_count = static_cast<std::size_t>(desc.halo) *
                       (shape.owned_nx + 2 * desc.halo) * shape.owned_nz;
  return shape;
}

void pack_face_x_faces(const state::FaceField<real>& field,
                       const domain::SubdomainDescriptor& desc,
                       ScalarHaloBuffers& buffers) {
  const auto shape = make_face_halo_shape(field, desc);
  const auto& storage = field.storage();
  std::size_t west_idx = 0;
  std::size_t east_idx = 0;

  for (int k = 0; k < shape.owned_nz; ++k) {
    for (int j = 0; j < shape.owned_ny; ++j) {
      for (int h = 0; h < desc.halo; ++h) {
        int west_i = h;
        int east_i = shape.owned_nx - desc.halo + h;
        if (field.orientation() == state::FaceOrientation::X) {
          west_i = 1 + h;
          east_i = desc.nx_local() - desc.halo + h;
        }
        buffers.west_send[west_idx++] = storage(west_i, j, k);
        buffers.east_send[east_idx++] = storage(east_i, j, k);
      }
    }
  }
}

void unpack_face_x_faces(state::FaceField<real>& field,
                         const domain::SubdomainDescriptor& desc,
                         const ScalarHaloBuffers& buffers) {
  const auto shape = make_face_halo_shape(field, desc);
  auto& storage = field.storage();
  std::size_t west_idx = 0;
  std::size_t east_idx = 0;
  const int east_halo_begin = shape.owned_nx;

  for (int k = 0; k < shape.owned_nz; ++k) {
    for (int j = 0; j < shape.owned_ny; ++j) {
      for (int h = 0; h < desc.halo; ++h) {
        storage(-desc.halo + h, j, k) = buffers.west_recv[west_idx++];
        storage(east_halo_begin + h, j, k) = buffers.east_recv[east_idx++];
      }
    }
  }
}

void pack_face_y_faces(const state::FaceField<real>& field,
                       const domain::SubdomainDescriptor& desc,
                       ScalarHaloBuffers& buffers) {
  const auto shape = make_face_halo_shape(field, desc);
  const auto& storage = field.storage();
  std::size_t south_idx = 0;
  std::size_t north_idx = 0;

  for (int k = 0; k < shape.owned_nz; ++k) {
    for (int h = 0; h < desc.halo; ++h) {
      for (int i = -desc.halo; i < shape.owned_nx + desc.halo; ++i) {
        int south_j = h;
        int north_j = shape.owned_ny - desc.halo + h;
        if (field.orientation() == state::FaceOrientation::Y) {
          south_j = 1 + h;
          north_j = desc.ny_local() - desc.halo + h;
        }
        buffers.south_send[south_idx++] = storage(i, south_j, k);
        buffers.north_send[north_idx++] = storage(i, north_j, k);
      }
    }
  }
}

void unpack_face_y_faces(state::FaceField<real>& field,
                         const domain::SubdomainDescriptor& desc,
                         const ScalarHaloBuffers& buffers) {
  const auto shape = make_face_halo_shape(field, desc);
  auto& storage = field.storage();
  std::size_t south_idx = 0;
  std::size_t north_idx = 0;
  const int north_halo_begin = shape.owned_ny;

  for (int k = 0; k < shape.owned_nz; ++k) {
    for (int h = 0; h < desc.halo; ++h) {
      for (int i = -desc.halo; i < shape.owned_nx + desc.halo; ++i) {
        storage(i, -desc.halo + h, k) = buffers.south_recv[south_idx++];
        storage(i, north_halo_begin + h, k) = buffers.north_recv[north_idx++];
      }
    }
  }
}

}  // namespace

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

void exchange_face_impl(
    const std::vector<state::FaceField<real>*>& fields,
    const std::vector<domain::SubdomainDescriptor>& layout) {
  gwm::require(fields.size() == layout.size(),
               "Field/layout size mismatch in exchange_face");
  if (fields.empty()) {
    return;
  }

  const auto reference_orientation = fields.front()->orientation();
  std::vector<ScalarHaloExchangePlan> plans;
  plans.reserve(layout.size());
  std::vector<ScalarHaloBuffers> buffers;
  buffers.reserve(layout.size());
  for (std::size_t n = 0; n < layout.size(); ++n) {
    gwm::require(fields[n] != nullptr, "Null field pointer in exchange_face");
    gwm::require(fields[n]->orientation() == reference_orientation,
                 "Face orientation mismatch in exchange_face");
    plans.push_back(make_scalar_halo_exchange_plan(layout, layout[n]));
    const auto shape = make_face_halo_shape(*fields[n], layout[n]);
    ScalarHaloExchangePlan plan = plans.back();
    plan.x_face_count = shape.x_face_count;
    plan.y_face_count = shape.y_face_count;
    buffers.push_back(HaloExchange::allocate_scalar_buffers(plan));
  }

  for (std::size_t n = 0; n < layout.size(); ++n) {
    pack_face_x_faces(*fields[n], layout[n], buffers[n]);
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
    unpack_face_x_faces(*fields[n], layout[n], buffers[n]);
  }

  for (std::size_t n = 0; n < layout.size(); ++n) {
    pack_face_y_faces(*fields[n], layout[n], buffers[n]);
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
    unpack_face_y_faces(*fields[n], layout[n], buffers[n]);
  }
}

void HaloExchange::exchange_face(
    std::vector<state::FaceField<real>>& fields,
    const std::vector<domain::SubdomainDescriptor>& layout) {
  std::vector<state::FaceField<real>*> field_ptrs;
  field_ptrs.reserve(fields.size());
  for (auto& field : fields) {
    field_ptrs.push_back(&field);
  }
  exchange_face_impl(field_ptrs, layout);
}

void HaloExchange::exchange_face(
    const std::vector<state::FaceField<real>*>& fields,
    const std::vector<domain::SubdomainDescriptor>& layout) {
  exchange_face_impl(fields, layout);
}

}  // namespace gwm::comm
