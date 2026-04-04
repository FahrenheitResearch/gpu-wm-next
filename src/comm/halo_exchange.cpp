#include "gwm/comm/halo_exchange.hpp"

#include <algorithm>

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

void HaloExchange::pack_face_x_faces(const state::FaceField<real>& field,
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

void HaloExchange::unpack_face_x_faces(state::FaceField<real>& field,
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

void HaloExchange::pack_face_y_faces(const state::FaceField<real>& field,
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

void HaloExchange::unpack_face_y_faces(state::FaceField<real>& field,
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

#if GWM_HAVE_MPI
MPI_Datatype mpi_real_datatype() { return MPI_FLOAT; }

void fill_recv_buffer(std::vector<real>& buffer) {
  std::fill(buffer.begin(), buffer.end(), 0.0f);
}

void mpi_sendrecv(const std::vector<real>& send_buffer, int dest_rank, int send_tag,
                  std::vector<real>& recv_buffer, int source_rank, int recv_tag,
                  const MpiCartesianContext& context) {
  fill_recv_buffer(recv_buffer);
  const int mpi_dest = dest_rank >= 0 ? dest_rank : MPI_PROC_NULL;
  const int mpi_source = source_rank >= 0 ? source_rank : MPI_PROC_NULL;
  MPI_Sendrecv(send_buffer.empty() ? nullptr : send_buffer.data(),
               static_cast<int>(send_buffer.size()), mpi_real_datatype(),
               mpi_dest, send_tag,
               recv_buffer.empty() ? nullptr : recv_buffer.data(),
               static_cast<int>(recv_buffer.size()), mpi_real_datatype(),
               mpi_source, recv_tag, context.cart_comm, MPI_STATUS_IGNORE);
}
#endif

void require_active_mpi_context(const MpiCartesianContext& context,
                                const domain::SubdomainDescriptor& desc) {
  gwm::require(context.active, "MPI halo exchange requires an active context");
  gwm::require(context.world_rank == desc.rank && context.coord_x == desc.coord_x &&
                   context.coord_y == desc.coord_y &&
                   context.ranks_x == desc.ranks_x &&
                   context.ranks_y == desc.ranks_y &&
                   context.periodic_x == desc.periodic_x &&
                   context.periodic_y == desc.periodic_y,
               "MPI cartesian context does not match the local descriptor");
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

ScalarHaloBuffers HaloExchange::allocate_face_buffers(
    const state::FaceField<real>& field,
    const domain::SubdomainDescriptor& desc) {
  auto plan = make_scalar_halo_exchange_plan(desc, CartesianNeighborRanks{});
  const auto shape = make_face_halo_shape(field, desc);
  plan.x_face_count = shape.x_face_count;
  plan.y_face_count = shape.y_face_count;
  return allocate_scalar_buffers(plan);
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

void HaloExchange::exchange_scalar(state::Field3D<real>& field,
                                   const domain::SubdomainDescriptor& desc,
                                   const MpiCartesianContext& context) {
  require_active_mpi_context(context, desc);
  const auto plan = make_scalar_halo_exchange_plan(desc, context.neighbors);
  auto buffers = allocate_scalar_buffers(plan);
  pack_scalar_x_faces(field, desc, buffers);
#if GWM_HAVE_MPI
  mpi_sendrecv(buffers.west_send, context.neighbors.west, 100, buffers.east_recv,
               context.neighbors.east, 100, context);
  mpi_sendrecv(buffers.east_send, context.neighbors.east, 101, buffers.west_recv,
               context.neighbors.west, 101, context);
#else
  (void)context;
  gwm::require(false, "MPI scalar halo exchange called without MPI support");
#endif
  unpack_scalar_x_faces(field, desc, buffers);

  pack_scalar_y_faces(field, desc, buffers);
#if GWM_HAVE_MPI
  mpi_sendrecv(buffers.south_send, context.neighbors.south, 110,
               buffers.north_recv, context.neighbors.north, 110, context);
  mpi_sendrecv(buffers.north_send, context.neighbors.north, 111,
               buffers.south_recv, context.neighbors.south, 111, context);
#endif
  unpack_scalar_y_faces(field, desc, buffers);
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

void HaloExchange::exchange_face(state::FaceField<real>& field,
                                 const domain::SubdomainDescriptor& desc,
                                 const MpiCartesianContext& context) {
  require_active_mpi_context(context, desc);

  ScalarHaloExchangePlan plan =
      make_scalar_halo_exchange_plan(desc, context.neighbors);
  const auto shape = make_face_halo_shape(field, desc);
  plan.x_face_count = shape.x_face_count;
  plan.y_face_count = shape.y_face_count;
  auto buffers = allocate_scalar_buffers(plan);

  pack_face_x_faces(field, desc, buffers);
#if GWM_HAVE_MPI
  mpi_sendrecv(buffers.west_send, context.neighbors.west, 120, buffers.east_recv,
               context.neighbors.east, 120, context);
  mpi_sendrecv(buffers.east_send, context.neighbors.east, 121, buffers.west_recv,
               context.neighbors.west, 121, context);
#else
  (void)context;
  gwm::require(false, "MPI face halo exchange called without MPI support");
#endif
  unpack_face_x_faces(field, desc, buffers);

  pack_face_y_faces(field, desc, buffers);
#if GWM_HAVE_MPI
  mpi_sendrecv(buffers.south_send, context.neighbors.south, 130,
               buffers.north_recv, context.neighbors.north, 130, context);
  mpi_sendrecv(buffers.north_send, context.neighbors.north, 131,
               buffers.south_recv, context.neighbors.south, 131, context);
#endif
  unpack_face_y_faces(field, desc, buffers);
}

void synchronize_face_owned_interfaces_impl(
    const std::vector<state::FaceField<real>*>& fields,
    const std::vector<domain::SubdomainDescriptor>& layout) {
  gwm::require(fields.size() == layout.size(),
               "Field/layout size mismatch in synchronize_owned_face_interfaces");
  if (fields.empty()) {
    return;
  }

  const auto orientation = fields.front()->orientation();
  if (orientation == state::FaceOrientation::Z) {
    return;
  }

  for (std::size_t n = 0; n < layout.size(); ++n) {
    gwm::require(fields[n] != nullptr,
                 "Null field pointer in synchronize_owned_face_interfaces");
    gwm::require(fields[n]->orientation() == orientation,
                 "Face orientation mismatch in synchronize_owned_face_interfaces");
  }

  for (std::size_t n = 0; n < layout.size(); ++n) {
    auto& field = *fields[n];
    auto& storage = field.storage();
    const auto& desc = layout[n];
    const auto neighbors = find_cartesian_neighbors(layout, desc);

    if (orientation == state::FaceOrientation::X && neighbors.east >= 0 &&
        (static_cast<int>(n) < neighbors.east ||
         neighbors.east == static_cast<int>(n))) {
      auto& peer_storage =
          fields[static_cast<std::size_t>(neighbors.east)]->storage();
      for (int k = 0; k < storage.nz(); ++k) {
        for (int j = 0; j < desc.ny_local(); ++j) {
          const real avg =
              0.5f * (storage(desc.nx_local(), j, k) + peer_storage(0, j, k));
          storage(desc.nx_local(), j, k) = avg;
          peer_storage(0, j, k) = avg;
        }
      }
    }

    if (orientation == state::FaceOrientation::Y && neighbors.north >= 0 &&
        (static_cast<int>(n) < neighbors.north ||
         neighbors.north == static_cast<int>(n))) {
      auto& peer_storage =
          fields[static_cast<std::size_t>(neighbors.north)]->storage();
      for (int k = 0; k < storage.nz(); ++k) {
        for (int i = 0; i < desc.nx_local(); ++i) {
          const real avg =
              0.5f * (storage(i, desc.ny_local(), k) + peer_storage(i, 0, k));
          storage(i, desc.ny_local(), k) = avg;
          peer_storage(i, 0, k) = avg;
        }
      }
    }
  }
}

void HaloExchange::synchronize_owned_face_interfaces(
    state::FaceField<real>& field, const domain::SubdomainDescriptor& desc,
    const MpiCartesianContext& context) {
  require_active_mpi_context(context, desc);

  if (field.orientation() == state::FaceOrientation::Z) {
    return;
  }

  ScalarHaloExchangePlan plan =
      make_scalar_halo_exchange_plan(desc, context.neighbors);
  const auto shape = make_face_halo_shape(field, desc);
  plan.x_face_count = shape.x_face_count;
  plan.y_face_count = shape.y_face_count;
  auto buffers = allocate_scalar_buffers(plan);
  auto& storage = field.storage();

#if GWM_HAVE_MPI
  if (field.orientation() == state::FaceOrientation::X) {
    const std::size_t interface_count =
        static_cast<std::size_t>(desc.ny_local()) * storage.nz();
    buffers.west_send.resize(interface_count);
    buffers.east_send.resize(interface_count);
    buffers.west_recv.resize(interface_count);
    buffers.east_recv.resize(interface_count);
    std::size_t west_idx = 0;
    std::size_t east_idx = 0;
    for (int k = 0; k < storage.nz(); ++k) {
      for (int j = 0; j < desc.ny_local(); ++j) {
        buffers.west_send[west_idx++] = storage(0, j, k);
        buffers.east_send[east_idx++] = storage(desc.nx_local(), j, k);
      }
    }
    mpi_sendrecv(buffers.west_send, context.neighbors.west, 140, buffers.west_recv,
                 context.neighbors.west, 141, context);
    mpi_sendrecv(buffers.east_send, context.neighbors.east, 141, buffers.east_recv,
                 context.neighbors.east, 140, context);
    west_idx = 0;
    east_idx = 0;
    if (context.neighbors.west >= 0) {
      for (int k = 0; k < storage.nz(); ++k) {
        for (int j = 0; j < desc.ny_local(); ++j) {
          storage(0, j, k) =
              0.5f * (storage(0, j, k) + buffers.west_recv[west_idx++]);
        }
      }
    }
    if (context.neighbors.east >= 0) {
      for (int k = 0; k < storage.nz(); ++k) {
        for (int j = 0; j < desc.ny_local(); ++j) {
          storage(desc.nx_local(), j, k) = 0.5f * (
              storage(desc.nx_local(), j, k) + buffers.east_recv[east_idx++]);
        }
      }
    }
    return;
  }

  const std::size_t interface_count =
      static_cast<std::size_t>(desc.nx_local()) * storage.nz();
  buffers.south_send.resize(interface_count);
  buffers.north_send.resize(interface_count);
  buffers.south_recv.resize(interface_count);
  buffers.north_recv.resize(interface_count);
  std::size_t south_idx = 0;
  std::size_t north_idx = 0;
  for (int k = 0; k < storage.nz(); ++k) {
    for (int i = 0; i < desc.nx_local(); ++i) {
      buffers.south_send[south_idx++] = storage(i, 0, k);
      buffers.north_send[north_idx++] = storage(i, desc.ny_local(), k);
    }
  }
  mpi_sendrecv(buffers.south_send, context.neighbors.south, 150,
               buffers.south_recv, context.neighbors.south, 151, context);
  mpi_sendrecv(buffers.north_send, context.neighbors.north, 151,
               buffers.north_recv, context.neighbors.north, 150, context);
  south_idx = 0;
  north_idx = 0;
  if (context.neighbors.south >= 0) {
    for (int k = 0; k < storage.nz(); ++k) {
      for (int i = 0; i < desc.nx_local(); ++i) {
        storage(i, 0, k) =
            0.5f * (storage(i, 0, k) + buffers.south_recv[south_idx++]);
      }
    }
  }
  if (context.neighbors.north >= 0) {
    for (int k = 0; k < storage.nz(); ++k) {
      for (int i = 0; i < desc.nx_local(); ++i) {
        storage(i, desc.ny_local(), k) = 0.5f * (
            storage(i, desc.ny_local(), k) + buffers.north_recv[north_idx++]);
      }
    }
  }
#else
  (void)context;
  gwm::require(false,
               "MPI face-interface synchronization called without MPI support");
#endif
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

void HaloExchange::synchronize_owned_face_interfaces(
    std::vector<state::FaceField<real>>& fields,
    const std::vector<domain::SubdomainDescriptor>& layout) {
  std::vector<state::FaceField<real>*> field_ptrs;
  field_ptrs.reserve(fields.size());
  for (auto& field : fields) {
    field_ptrs.push_back(&field);
  }
  synchronize_face_owned_interfaces_impl(field_ptrs, layout);
}

void HaloExchange::synchronize_owned_face_interfaces(
    const std::vector<state::FaceField<real>*>& fields,
    const std::vector<domain::SubdomainDescriptor>& layout) {
  synchronize_face_owned_interfaces_impl(fields, layout);
}

}  // namespace gwm::comm
