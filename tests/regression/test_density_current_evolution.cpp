#include <cmath>
#include <string>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/dry_diagnostics.hpp"
#include "gwm/dycore/idealized_cases.hpp"

#include "test_assert.hpp"

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig cfg{};
    cfg.nx = 32;
    cfg.ny = 32;
    cfg.nz = 16;
    cfg.halo = 1;
    cfg.ranks_x = 2;
    cfg.ranks_y = 2;

    const auto domain = domain::build_rectilinear_domain(cfg);

    dycore::ThermoBubbleConfig bubble_cfg{};
    bubble_cfg.rho_background = 1.0f;
    bubble_cfg.theta_background = 300.0f;
    bubble_cfg.theta_perturbation = 3.0f;
    bubble_cfg.center_x_fraction = 0.5f;
    bubble_cfg.center_y_fraction = 0.5f;
    bubble_cfg.center_z_fraction = 0.30f;
    bubble_cfg.radius_x_fraction = 0.12f;
    bubble_cfg.radius_y_fraction = 0.12f;
    bubble_cfg.radius_z_fraction = 0.14f;

    auto states = dycore::make_density_current_state(
        domain.layout, domain.metrics, bubble_cfg, "density_current");
    const auto initial = dycore::summarize_dry_states(states);

    dycore::DryStepperConfig step_cfg{};
    step_cfg.dt = 2.0f;
    step_cfg.fast_substeps = 3;
    dycore::NullBoundaryUpdater boundary;
    dycore::LocalSplitExplicitFastMode fast_modes;

    for (int step = 0; step < 5; ++step) {
      dycore::advance_dry_state_ssprk3(states, domain.layout, domain.metrics,
                                       step_cfg, boundary, fast_modes);
    }

    const auto final = dycore::summarize_dry_states(states);
    TEST_NEAR(initial.total_dry_mass, final.total_dry_mass, 2.0e-2);
    TEST_CHECK(final.u_face.max > 1.0e-3);
    TEST_CHECK(final.u_face.min < -1.0e-3);

    const auto rho = comm::VirtualRankLayout::gather_scalar(
        [&states]() {
          std::vector<state::Field3D<real>> fields;
          for (const auto& s : states) {
            fields.push_back(s.rho_d.clone_empty_like("_rho_g"));
            fields.back().copy_all_from(s.rho_d);
          }
          return fields;
        }(),
        domain.layout, "rho_density_current");

    const int ic = rho.nx() / 2;
    const int jc = rho.ny() / 2;
    const int kc = rho.nz() / 3;

    for (int offset = 1; offset <= 4; ++offset) {
      const double left = static_cast<double>(rho(ic - offset, jc, kc));
      const double right = static_cast<double>(rho(ic + offset, jc, kc));
      if (std::fabs(left - right) > 7.0e-3) {
        test_fail("rho symmetry offset " + std::to_string(offset) + " diff=" +
                  std::to_string(std::fabs(left - right)));
      }
    }

    for (int k = 0; k < rho.nz(); ++k) {
      for (int j = 0; j < rho.ny(); ++j) {
        for (int i = 0; i < rho.nx(); ++i) {
          TEST_CHECK(std::isfinite(static_cast<double>(rho(i, j, k))));
        }
      }
    }

    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
