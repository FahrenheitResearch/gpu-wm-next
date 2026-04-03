#include "gwm/comm/virtual_rank_layout.hpp"

#include <string>

namespace gwm::comm {

std::vector<domain::SubdomainDescriptor> VirtualRankLayout::build(
    int global_nx, int global_ny, int nz, int halo, int ranks_x, int ranks_y,
    bool periodic_x, bool periodic_y) {
  std::vector<domain::SubdomainDescriptor> layout;
  layout.reserve(static_cast<std::size_t>(ranks_x * ranks_y));

  const int base_nx = global_nx / ranks_x;
  const int rem_nx = global_nx % ranks_x;
  const int base_ny = global_ny / ranks_y;
  const int rem_ny = global_ny % ranks_y;

  int j_begin = 0;
  for (int ry = 0; ry < ranks_y; ++ry) {
    const int ny_local = base_ny + (ry < rem_ny ? 1 : 0);
    int i_begin = 0;
    for (int rx = 0; rx < ranks_x; ++rx) {
      const int nx_local = base_nx + (rx < rem_nx ? 1 : 0);
      domain::SubdomainDescriptor desc{};
      desc.rank = static_cast<int>(layout.size());
      desc.coord_x = rx;
      desc.coord_y = ry;
      desc.ranks_x = ranks_x;
      desc.ranks_y = ranks_y;
      desc.global_nx = global_nx;
      desc.global_ny = global_ny;
      desc.nz = nz;
      desc.halo = halo;
      desc.i_begin = i_begin;
      desc.i_end = i_begin + nx_local;
      desc.j_begin = j_begin;
      desc.j_end = j_begin + ny_local;
      desc.periodic_x = periodic_x;
      desc.periodic_y = periodic_y;
      layout.push_back(desc);
      i_begin += nx_local;
    }
    j_begin += ny_local;
  }

  return layout;
}

std::vector<state::Field3D<real>> VirtualRankLayout::scatter_scalar(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const std::function<real(int, int, int)>& initializer,
    const std::string& label_prefix) {
  std::vector<state::Field3D<real>> fields;
  fields.reserve(layout.size());
  for (const auto& desc : layout) {
    state::Field3D<real> field(desc.nx_local(), desc.ny_local(), desc.nz,
                               desc.halo,
                               label_prefix + "_rank_" +
                                   std::to_string(desc.rank));
    field.fill(0.0f);
    for (int k = 0; k < desc.nz; ++k) {
      for (int j = 0; j < desc.ny_local(); ++j) {
        for (int i = 0; i < desc.nx_local(); ++i) {
          field(i, j, k) = initializer(desc.i_begin + i, desc.j_begin + j, k);
        }
      }
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

state::Field3D<real> VirtualRankLayout::gather_scalar(
    const std::vector<state::Field3D<real>>& fields,
    const std::vector<domain::SubdomainDescriptor>& layout, const std::string& label) {
  gwm::require(!layout.empty(), "Layout must not be empty in gather_scalar");
  state::Field3D<real> global(layout.front().global_nx, layout.front().global_ny,
                              layout.front().nz, 0, label);
  global.fill(0.0f);

  for (std::size_t n = 0; n < layout.size(); ++n) {
    const auto& desc = layout[n];
    const auto& local = fields[n];
    for (int k = 0; k < desc.nz; ++k) {
      for (int j = 0; j < desc.ny_local(); ++j) {
        for (int i = 0; i < desc.nx_local(); ++i) {
          global(desc.i_begin + i, desc.j_begin + j, k) = local(i, j, k);
        }
      }
    }
  }

  return global;
}

state::FaceField<real> VirtualRankLayout::gather_face(
    const std::vector<state::FaceField<real>>& fields,
    const std::vector<domain::SubdomainDescriptor>& layout,
    state::FaceOrientation orientation, const std::string& label) {
  gwm::require(!layout.empty(), "Layout must not be empty in gather_face");

  const int global_nx =
      layout.front().global_nx +
      (orientation == state::FaceOrientation::X ? 1 : 0);
  const int global_ny =
      layout.front().global_ny +
      (orientation == state::FaceOrientation::Y ? 1 : 0);
  const int global_nz =
      layout.front().nz + (orientation == state::FaceOrientation::Z ? 1 : 0);

  state::FaceField<real> global(layout.front().global_nx, layout.front().global_ny,
                                layout.front().nz, 0, orientation, label);
  global.storage().fill(0.0f);

  for (std::size_t n = 0; n < layout.size(); ++n) {
    const auto& desc = layout[n];
    const auto& local = fields[n].storage();
    gwm::require(fields[n].orientation() == orientation,
                 "Face orientation mismatch in gather_face");
    for (int k = 0; k < local.nz(); ++k) {
      for (int j = 0; j < local.ny(); ++j) {
        for (int i = 0; i < local.nx(); ++i) {
          const int gi = desc.i_begin + i;
          const int gj = desc.j_begin + j;
          gwm::require(gi >= 0 && gi < global_nx && gj >= 0 && gj < global_ny &&
                           k >= 0 && k < global_nz,
                       "Out-of-range face gather index");
          global.storage()(gi, gj, k) = local(i, j, k);
        }
      }
    }
  }

  return global;
}

}  // namespace gwm::comm
