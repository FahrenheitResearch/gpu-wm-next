#include "gwm/dycore/boundary_conditions.hpp"

#include <string>

namespace gwm::dycore {

namespace {

void copy_west_cell_strip(state::Field3D<real>& target,
                          const state::Field3D<real>& source) {
  for (int k = 0; k < target.nz(); ++k) {
    for (int j = 0; j < target.ny(); ++j) {
      target(0, j, k) = source(0, j, k);
    }
  }
}

void copy_east_cell_strip(state::Field3D<real>& target,
                          const state::Field3D<real>& source) {
  const int i = target.nx() - 1;
  for (int k = 0; k < target.nz(); ++k) {
    for (int j = 0; j < target.ny(); ++j) {
      target(i, j, k) = source(i, j, k);
    }
  }
}

void copy_south_cell_strip(state::Field3D<real>& target,
                           const state::Field3D<real>& source) {
  for (int k = 0; k < target.nz(); ++k) {
    for (int i = 0; i < target.nx(); ++i) {
      target(i, 0, k) = source(i, 0, k);
    }
  }
}

void copy_north_cell_strip(state::Field3D<real>& target,
                           const state::Field3D<real>& source) {
  const int j = target.ny() - 1;
  for (int k = 0; k < target.nz(); ++k) {
    for (int i = 0; i < target.nx(); ++i) {
      target(i, j, k) = source(i, j, k);
    }
  }
}

void copy_west_face_strip(state::FaceField<real>& target,
                          const state::FaceField<real>& source) {
  auto& target_storage = target.storage();
  const auto& source_storage = source.storage();
  for (int k = 0; k < target_storage.nz(); ++k) {
    for (int j = 0; j < target_storage.ny(); ++j) {
      target_storage(0, j, k) = source_storage(0, j, k);
    }
  }
}

void copy_east_face_strip(state::FaceField<real>& target,
                          const state::FaceField<real>& source) {
  auto& target_storage = target.storage();
  const auto& source_storage = source.storage();
  const int i = target_storage.nx() - 1;
  for (int k = 0; k < target_storage.nz(); ++k) {
    for (int j = 0; j < target_storage.ny(); ++j) {
      target_storage(i, j, k) = source_storage(i, j, k);
    }
  }
}

void copy_south_face_strip(state::FaceField<real>& target,
                           const state::FaceField<real>& source) {
  auto& target_storage = target.storage();
  const auto& source_storage = source.storage();
  for (int k = 0; k < target_storage.nz(); ++k) {
    for (int i = 0; i < target_storage.nx(); ++i) {
      target_storage(i, 0, k) = source_storage(i, 0, k);
    }
  }
}

void copy_north_face_strip(state::FaceField<real>& target,
                           const state::FaceField<real>& source) {
  auto& target_storage = target.storage();
  const auto& source_storage = source.storage();
  const int j = target_storage.ny() - 1;
  for (int k = 0; k < target_storage.nz(); ++k) {
    for (int i = 0; i < target_storage.nx(); ++i) {
      target_storage(i, j, k) = source_storage(i, j, k);
    }
  }
}

}  // namespace

bool touches_west_boundary(const domain::SubdomainDescriptor& desc) {
  return !desc.periodic_x && desc.coord_x == 0;
}

bool touches_east_boundary(const domain::SubdomainDescriptor& desc) {
  return !desc.periodic_x && desc.coord_x == desc.ranks_x - 1;
}

bool touches_south_boundary(const domain::SubdomainDescriptor& desc) {
  return !desc.periodic_y && desc.coord_y == 0;
}

bool touches_north_boundary(const domain::SubdomainDescriptor& desc) {
  return !desc.periodic_y && desc.coord_y == desc.ranks_y - 1;
}

void apply_reference_boundaries(
    std::vector<DryState>& states, const std::vector<DryState>& boundary_states,
    const std::vector<domain::SubdomainDescriptor>& layout) {
  gwm::require(states.size() == boundary_states.size(),
               "State/boundary-state size mismatch in apply_reference_boundaries");
  gwm::require(states.size() == layout.size(),
               "State/layout size mismatch in apply_reference_boundaries");

  for (std::size_t idx = 0; idx < states.size(); ++idx) {
    auto& state = states[idx];
    const auto& boundary = boundary_states[idx];
    const auto& desc = layout[idx];

    gwm::require(state.rho_d.nx() == boundary.rho_d.nx() &&
                     state.rho_d.ny() == boundary.rho_d.ny() &&
                     state.rho_d.nz() == boundary.rho_d.nz(),
                 "Boundary reference state dimensions must match target state");

    if (touches_west_boundary(desc)) {
      copy_west_cell_strip(state.rho_d, boundary.rho_d);
      copy_west_cell_strip(state.rho_theta_m, boundary.rho_theta_m);
      copy_west_face_strip(state.mom_u, boundary.mom_u);
      copy_west_face_strip(state.mom_v, boundary.mom_v);
      copy_west_face_strip(state.mom_w, boundary.mom_w);
    }

    if (touches_east_boundary(desc)) {
      copy_east_cell_strip(state.rho_d, boundary.rho_d);
      copy_east_cell_strip(state.rho_theta_m, boundary.rho_theta_m);
      copy_east_face_strip(state.mom_u, boundary.mom_u);
      copy_east_face_strip(state.mom_v, boundary.mom_v);
      copy_east_face_strip(state.mom_w, boundary.mom_w);
    }

    if (touches_south_boundary(desc)) {
      copy_south_cell_strip(state.rho_d, boundary.rho_d);
      copy_south_cell_strip(state.rho_theta_m, boundary.rho_theta_m);
      copy_south_face_strip(state.mom_u, boundary.mom_u);
      copy_south_face_strip(state.mom_v, boundary.mom_v);
      copy_south_face_strip(state.mom_w, boundary.mom_w);
    }

    if (touches_north_boundary(desc)) {
      copy_north_cell_strip(state.rho_d, boundary.rho_d);
      copy_north_cell_strip(state.rho_theta_m, boundary.rho_theta_m);
      copy_north_face_strip(state.mom_u, boundary.mom_u);
      copy_north_face_strip(state.mom_v, boundary.mom_v);
      copy_north_face_strip(state.mom_w, boundary.mom_w);
    }
  }
}

void apply_reference_boundaries(
    std::vector<state::TracerState>& states,
    const std::vector<state::TracerState>& boundary_states,
    const std::vector<domain::SubdomainDescriptor>& layout) {
  gwm::require(states.size() == boundary_states.size(),
               "Tracer-state/boundary-state size mismatch in "
               "apply_reference_boundaries");
  gwm::require(states.size() == layout.size(),
               "Tracer-state/layout size mismatch in "
               "apply_reference_boundaries");

  for (std::size_t idx = 0; idx < states.size(); ++idx) {
    auto& state = states[idx];
    const auto& boundary = boundary_states[idx];
    const auto& desc = layout[idx];

    gwm::require(state.size() == boundary.size(),
                 "Tracer-state count mismatch in apply_reference_boundaries");
    for (std::size_t tracer = 0; tracer < state.size(); ++tracer) {
      auto& target_field = state.mass(static_cast<int>(tracer));
      const auto& source_field = boundary.mass(static_cast<int>(tracer));
      gwm::require(target_field.nx() == source_field.nx() &&
                       target_field.ny() == source_field.ny() &&
                       target_field.nz() == source_field.nz(),
                   "Tracer boundary reference dimensions must match target "
                   "state");

      if (touches_west_boundary(desc)) {
        copy_west_cell_strip(target_field, source_field);
      }

      if (touches_east_boundary(desc)) {
        copy_east_cell_strip(target_field, source_field);
      }

      if (touches_south_boundary(desc)) {
        copy_south_cell_strip(target_field, source_field);
      }

      if (touches_north_boundary(desc)) {
        copy_north_cell_strip(target_field, source_field);
      }
    }
  }
}

}  // namespace gwm::dycore
