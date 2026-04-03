#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"

#include "test_assert.hpp"

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig domain_cfg{};
    domain_cfg.nx = 20;
    domain_cfg.ny = 16;
    domain_cfg.nz = 6;
    domain_cfg.halo = 1;
    domain_cfg.ranks_x = 2;
    domain_cfg.ranks_y = 2;

    const auto domain = domain::build_rectilinear_domain(domain_cfg);
    auto states = dycore::make_constant_dry_state(
        domain.layout, 1.1f, 301.0f, 10.0f, -4.0f, 0.0f, "const");

    dycore::DryStepperConfig step_cfg{};
    step_cfg.dt = 2.0f;
    step_cfg.fast_substeps = 2;
    dycore::NullBoundaryUpdater boundary;
    dycore::LocalSplitExplicitFastMode fast_modes;
    dycore::advance_dry_state_ssprk3(states, domain.layout, domain.metrics,
                                     step_cfg, boundary, fast_modes);

    for (const auto& state : states) {
      for (int k = 0; k < state.rho_d.nz(); ++k) {
        for (int j = 0; j < state.rho_d.ny(); ++j) {
          for (int i = 0; i < state.rho_d.nx(); ++i) {
            TEST_NEAR(state.rho_d(i, j, k), 1.1f, 1.0e-5f);
            TEST_NEAR(state.rho_theta_m(i, j, k), 1.1f * 301.0f, 1.0e-3f);
          }
        }
      }
    }
    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
