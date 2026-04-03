#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/dry_diagnostics.hpp"

#include "test_assert.hpp"

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig cfg_domain{};
    cfg_domain.nx = 18;
    cfg_domain.ny = 18;
    cfg_domain.nz = 12;
    cfg_domain.halo = 1;
    cfg_domain.ranks_x = 2;
    cfg_domain.ranks_y = 2;
    cfg_domain.dx = 3000.0;
    cfg_domain.dy = 3000.0;
    cfg_domain.z_top = 18000.0;

    const auto domain = domain::build_rectilinear_domain(cfg_domain);
    auto states = dycore::make_hydrostatic_rest_state(
        domain.layout, domain.metrics, 1.22f, 300.0f, 8200.0f, "hydro");

    double mass_before = 0.0;
    double theta_before = 0.0;
    for (const auto& state : states) {
      mass_before += state.total_dry_mass();
      theta_before += state.total_rho_theta_m();
    }

    dycore::DryStepperConfig step_cfg{};
    step_cfg.dt = 1.0f;
    step_cfg.fast_substeps = 3;
    dycore::NullBoundaryUpdater boundary;
    dycore::LocalSplitExplicitFastMode fast_modes;

    dycore::advance_dry_state_ssprk3(states, domain.layout, domain.metrics,
                                     step_cfg, boundary, fast_modes);

    double mass_after = 0.0;
    double theta_after = 0.0;
    for (const auto& state : states) {
      mass_after += state.total_dry_mass();
      theta_after += state.total_rho_theta_m();
    }
    const auto summary = dycore::summarize_dry_states(states);

    TEST_NEAR(mass_before, mass_after, 1.0e-3);
    TEST_NEAR(theta_before, theta_after, 3.0e-2);
    TEST_NEAR(summary.total_horizontal_momentum_x, 0.0, 1.0e-4);
    TEST_NEAR(summary.total_horizontal_momentum_y, 0.0, 1.0e-4);
    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
