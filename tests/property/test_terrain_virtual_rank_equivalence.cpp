#include <algorithm>
#include <cmath>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/idealized_cases.hpp"

#include "test_assert.hpp"

namespace {

void compare_scalar_fields(const gwm::state::Field3D<gwm::real>& a,
                           const gwm::state::Field3D<gwm::real>& b,
                           double tol, const char* label) {
  TEST_CHECK(a.nx() == b.nx());
  TEST_CHECK(a.ny() == b.ny());
  TEST_CHECK(a.nz() == b.nz());
  double max_diff = 0.0;
  for (int k = 0; k < a.nz(); ++k) {
    for (int j = 0; j < a.ny(); ++j) {
      for (int i = 0; i < a.nx(); ++i) {
        const double diff = std::fabs(static_cast<double>(a(i, j, k)) -
                                      static_cast<double>(b(i, j, k)));
        max_diff = std::max(max_diff, diff);
      }
    }
  }
  if (max_diff > tol) {
    test_fail(std::string(label) + " max diff " + std::to_string(max_diff));
  }
}

}  // namespace

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig serial_cfg{};
    serial_cfg.nx = 64;
    serial_cfg.ny = 24;
    serial_cfg.nz = 20;
    serial_cfg.halo = 1;
    serial_cfg.ranks_x = 1;
    serial_cfg.ranks_y = 1;
    serial_cfg.dx = 1000.0;
    serial_cfg.dy = 1000.0;
    serial_cfg.z_top = 18000.0;
    serial_cfg.terrain_kind = domain::TerrainProfileKind::CosineMountain;
    serial_cfg.mountain_height = 1000.0;
    serial_cfg.mountain_half_width_x = 0.12;
    serial_cfg.mountain_half_width_y = 0.40;

    domain::RectilinearDomainConfig split_cfg = serial_cfg;
    split_cfg.ranks_x = 2;
    split_cfg.ranks_y = 2;

    const auto serial_domain = domain::build_rectilinear_domain(serial_cfg);
    const auto split_domain = domain::build_rectilinear_domain(split_cfg);

    dycore::MountainWaveBackgroundConfig bg{};
    bg.rho_surface = 1.18f;
    bg.theta_ref = 300.0f;
    bg.density_scale_height = 8200.0f;
    bg.u_background = 15.0f;

    auto serial_state = dycore::make_mountain_wave_background_state(
        serial_domain.layout, serial_domain.metrics, serial_domain, bg, "serial_terrain");
    auto split_state = dycore::make_mountain_wave_background_state(
        split_domain.layout, split_domain.metrics, split_domain, bg, "split_terrain");

    dycore::DryStepperConfig cfg{};
    cfg.dt = 1.0f;
    cfg.fast_substeps = 4;
    dycore::NullBoundaryUpdater boundary;
    dycore::LocalSplitExplicitFastMode fast_modes;

    for (int step = 0; step < 3; ++step) {
      dycore::advance_dry_state_ssprk3(serial_state, serial_domain.layout,
                                       serial_domain.metrics, cfg, boundary,
                                       fast_modes);
      dycore::advance_dry_state_ssprk3(split_state, split_domain.layout,
                                       split_domain.metrics, cfg, boundary,
                                       fast_modes);
    }

    const auto serial_rho = comm::VirtualRankLayout::gather_scalar(
        [&serial_state]() {
          std::vector<state::Field3D<real>> fields;
          for (const auto& s : serial_state) {
            fields.push_back(s.rho_d.clone_empty_like("_g"));
            fields.back().copy_all_from(s.rho_d);
          }
          return fields;
        }(),
        serial_domain.layout, "serial_terrain_rho");
    const auto split_rho = comm::VirtualRankLayout::gather_scalar(
        [&split_state]() {
          std::vector<state::Field3D<real>> fields;
          for (const auto& s : split_state) {
            fields.push_back(s.rho_d.clone_empty_like("_g"));
            fields.back().copy_all_from(s.rho_d);
          }
          return fields;
        }(),
        split_domain.layout, "split_terrain_rho");

    const auto serial_theta = comm::VirtualRankLayout::gather_scalar(
        [&serial_state]() {
          std::vector<state::Field3D<real>> fields;
          for (const auto& s : serial_state) {
            fields.push_back(s.rho_theta_m.clone_empty_like("_g"));
            fields.back().copy_all_from(s.rho_theta_m);
          }
          return fields;
        }(),
        serial_domain.layout, "serial_terrain_theta");
    const auto split_theta = comm::VirtualRankLayout::gather_scalar(
        [&split_state]() {
          std::vector<state::Field3D<real>> fields;
          for (const auto& s : split_state) {
            fields.push_back(s.rho_theta_m.clone_empty_like("_g"));
            fields.back().copy_all_from(s.rho_theta_m);
          }
          return fields;
        }(),
        split_domain.layout, "split_terrain_theta");

    const auto serial_u = comm::VirtualRankLayout::gather_face(
        [&serial_state]() {
          std::vector<state::FaceField<real>> fields;
          for (const auto& s : serial_state) {
            fields.push_back(state::FaceField<real>(s.rho_d.nx(), s.rho_d.ny(),
                                                    s.rho_d.nz(), s.rho_d.halo(),
                                                    state::FaceOrientation::X, "u"));
            fields.back().storage().copy_all_from(s.mom_u.storage());
          }
          return fields;
        }(),
        serial_domain.layout, state::FaceOrientation::X, "serial_terrain_u");
    const auto split_u = comm::VirtualRankLayout::gather_face(
        [&split_state]() {
          std::vector<state::FaceField<real>> fields;
          for (const auto& s : split_state) {
            fields.push_back(state::FaceField<real>(s.rho_d.nx(), s.rho_d.ny(),
                                                    s.rho_d.nz(), s.rho_d.halo(),
                                                    state::FaceOrientation::X, "u"));
            fields.back().storage().copy_all_from(s.mom_u.storage());
          }
          return fields;
        }(),
        split_domain.layout, state::FaceOrientation::X, "split_terrain_u");

    const auto serial_v = comm::VirtualRankLayout::gather_face(
        [&serial_state]() {
          std::vector<state::FaceField<real>> fields;
          for (const auto& s : serial_state) {
            fields.push_back(state::FaceField<real>(s.rho_d.nx(), s.rho_d.ny(),
                                                    s.rho_d.nz(), s.rho_d.halo(),
                                                    state::FaceOrientation::Y, "v"));
            fields.back().storage().copy_all_from(s.mom_v.storage());
          }
          return fields;
        }(),
        serial_domain.layout, state::FaceOrientation::Y, "serial_terrain_v");
    const auto split_v = comm::VirtualRankLayout::gather_face(
        [&split_state]() {
          std::vector<state::FaceField<real>> fields;
          for (const auto& s : split_state) {
            fields.push_back(state::FaceField<real>(s.rho_d.nx(), s.rho_d.ny(),
                                                    s.rho_d.nz(), s.rho_d.halo(),
                                                    state::FaceOrientation::Y, "v"));
            fields.back().storage().copy_all_from(s.mom_v.storage());
          }
          return fields;
        }(),
        split_domain.layout, state::FaceOrientation::Y, "split_terrain_v");

    const auto serial_w = comm::VirtualRankLayout::gather_face(
        [&serial_state]() {
          std::vector<state::FaceField<real>> fields;
          for (const auto& s : serial_state) {
            fields.push_back(state::FaceField<real>(s.rho_d.nx(), s.rho_d.ny(),
                                                    s.rho_d.nz(), s.rho_d.halo(),
                                                    state::FaceOrientation::Z, "w"));
            fields.back().storage().copy_all_from(s.mom_w.storage());
          }
          return fields;
        }(),
        serial_domain.layout, state::FaceOrientation::Z, "serial_terrain_w");
    const auto split_w = comm::VirtualRankLayout::gather_face(
        [&split_state]() {
          std::vector<state::FaceField<real>> fields;
          for (const auto& s : split_state) {
            fields.push_back(state::FaceField<real>(s.rho_d.nx(), s.rho_d.ny(),
                                                    s.rho_d.nz(), s.rho_d.halo(),
                                                    state::FaceOrientation::Z, "w"));
            fields.back().storage().copy_all_from(s.mom_w.storage());
          }
          return fields;
        }(),
        split_domain.layout, state::FaceOrientation::Z, "split_terrain_w");

    compare_scalar_fields(serial_rho, split_rho, 2.0e-4, "terrain_rho");
    compare_scalar_fields(serial_theta, split_theta, 1.0e-3, "terrain_theta");
    compare_scalar_fields(serial_u.storage(), split_u.storage(), 3.0e-4, "terrain_u");
    compare_scalar_fields(serial_v.storage(), split_v.storage(), 3.0e-4, "terrain_v");
    compare_scalar_fields(serial_w.storage(), split_w.storage(), 5.0e-4, "terrain_w");
    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
