#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/idealized_cases.hpp"

#include "test_assert.hpp"

namespace {

void compare_scalar_fields(const gwm::state::Field3D<gwm::real>& a,
                           const gwm::state::Field3D<gwm::real>& b,
                           double tol) {
  TEST_CHECK(a.nx() == b.nx());
  TEST_CHECK(a.ny() == b.ny());
  TEST_CHECK(a.nz() == b.nz());
  for (int k = 0; k < a.nz(); ++k) {
    for (int j = 0; j < a.ny(); ++j) {
      for (int i = 0; i < a.nx(); ++i) {
        TEST_NEAR(a(i, j, k), b(i, j, k), tol);
      }
    }
  }
}

}  // namespace

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig serial_cfg{};
    serial_cfg.nx = 24;
    serial_cfg.ny = 20;
    serial_cfg.nz = 8;
    serial_cfg.halo = 1;
    serial_cfg.ranks_x = 1;
    serial_cfg.ranks_y = 1;

    domain::RectilinearDomainConfig split_cfg = serial_cfg;
    split_cfg.ranks_x = 2;
    split_cfg.ranks_y = 2;

    const auto serial_domain = domain::build_rectilinear_domain(serial_cfg);
    const auto split_domain = domain::build_rectilinear_domain(split_cfg);

    dycore::ThermoBubbleConfig bubble_cfg{};
    bubble_cfg.rho_background = 1.0f;
    bubble_cfg.theta_background = 300.0f;
    bubble_cfg.theta_perturbation = 2.5f;
    bubble_cfg.center_x_fraction = 0.5f;
    bubble_cfg.center_y_fraction = 0.5f;
    bubble_cfg.center_z_fraction = 0.35f;
    bubble_cfg.radius_x_fraction = 0.12f;
    bubble_cfg.radius_y_fraction = 0.12f;
    bubble_cfg.radius_z_fraction = 0.16f;

    auto serial_state = dycore::make_warm_bubble_state(
        serial_domain.layout, serial_domain.metrics, bubble_cfg, "serial_fast");
    auto split_state = dycore::make_warm_bubble_state(
        split_domain.layout, split_domain.metrics, bubble_cfg, "split_fast");

    dycore::DryStepperConfig cfg{};
    cfg.dt = 1.0f;
    cfg.fast_substeps = 4;
    dycore::LocalSplitExplicitFastMode fast_modes;

    fast_modes.apply(serial_state, serial_domain.layout, serial_domain.metrics,
                     cfg);
    fast_modes.apply(split_state, split_domain.layout, split_domain.metrics,
                     cfg);

    const auto serial_rho =
        comm::VirtualRankLayout::gather_scalar(
            [&serial_state]() {
              std::vector<state::Field3D<real>> fields;
              for (const auto& s : serial_state) {
                fields.push_back(s.rho_d.clone_empty_like("_g"));
                fields.back().copy_all_from(s.rho_d);
              }
              return fields;
            }(),
            serial_domain.layout, "serial_fast_rho");
    const auto split_rho =
        comm::VirtualRankLayout::gather_scalar(
            [&split_state]() {
              std::vector<state::Field3D<real>> fields;
              for (const auto& s : split_state) {
                fields.push_back(s.rho_d.clone_empty_like("_g"));
                fields.back().copy_all_from(s.rho_d);
              }
              return fields;
            }(),
            split_domain.layout, "split_fast_rho");

    const auto serial_u = comm::VirtualRankLayout::gather_face(
        [&serial_state]() {
          std::vector<state::FaceField<real>> fields;
          for (const auto& s : serial_state) {
            fields.push_back(state::FaceField<real>(s.rho_d.nx(), s.rho_d.ny(),
                                                    s.rho_d.nz(), s.rho_d.halo(),
                                                    state::FaceOrientation::X,
                                                    "u"));
            fields.back().storage().copy_all_from(s.mom_u.storage());
          }
          return fields;
        }(),
        serial_domain.layout, state::FaceOrientation::X, "serial_fast_u");
    const auto split_u = comm::VirtualRankLayout::gather_face(
        [&split_state]() {
          std::vector<state::FaceField<real>> fields;
          for (const auto& s : split_state) {
            fields.push_back(state::FaceField<real>(s.rho_d.nx(), s.rho_d.ny(),
                                                    s.rho_d.nz(), s.rho_d.halo(),
                                                    state::FaceOrientation::X,
                                                    "u"));
            fields.back().storage().copy_all_from(s.mom_u.storage());
          }
          return fields;
        }(),
        split_domain.layout, state::FaceOrientation::X, "split_fast_u");

    const auto serial_v = comm::VirtualRankLayout::gather_face(
        [&serial_state]() {
          std::vector<state::FaceField<real>> fields;
          for (const auto& s : serial_state) {
            fields.push_back(state::FaceField<real>(s.rho_d.nx(), s.rho_d.ny(),
                                                    s.rho_d.nz(), s.rho_d.halo(),
                                                    state::FaceOrientation::Y,
                                                    "v"));
            fields.back().storage().copy_all_from(s.mom_v.storage());
          }
          return fields;
        }(),
        serial_domain.layout, state::FaceOrientation::Y, "serial_fast_v");
    const auto split_v = comm::VirtualRankLayout::gather_face(
        [&split_state]() {
          std::vector<state::FaceField<real>> fields;
          for (const auto& s : split_state) {
            fields.push_back(state::FaceField<real>(s.rho_d.nx(), s.rho_d.ny(),
                                                    s.rho_d.nz(), s.rho_d.halo(),
                                                    state::FaceOrientation::Y,
                                                    "v"));
            fields.back().storage().copy_all_from(s.mom_v.storage());
          }
          return fields;
        }(),
        split_domain.layout, state::FaceOrientation::Y, "split_fast_v");

    compare_scalar_fields(serial_rho, split_rho, 1.0e-4);
    compare_scalar_fields(serial_u.storage(), split_u.storage(), 1.0e-4);
    compare_scalar_fields(serial_v.storage(), split_v.storage(), 1.0e-4);
    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
