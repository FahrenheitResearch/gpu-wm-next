#include <cmath>

#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/dry_diagnostics.hpp"

#include "test_assert.hpp"

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig cfg{};
    cfg.nx = 33;
    cfg.ny = 33;
    cfg.nz = 9;
    cfg.halo = 1;
    cfg.ranks_x = 1;
    cfg.ranks_y = 1;

    const auto domain = domain::build_rectilinear_domain(cfg);
    auto states = dycore::make_constant_dry_state(
        domain.layout, 1.0f, 300.0f, 0.0f, 0.0f, 0.0f, "acoustic");

    auto& state = states.front();
    const int ic = cfg.nx / 2;
    const int jc = cfg.ny / 2;
    const int kc = cfg.nz / 2;
    for (int k = 0; k < state.rho_d.nz(); ++k) {
      for (int j = 0; j < state.rho_d.ny(); ++j) {
        for (int i = 0; i < state.rho_d.nx(); ++i) {
          const real dx = static_cast<real>(i - ic);
          const real dy = static_cast<real>(j - jc);
          const real dz = static_cast<real>(k - kc);
          const real r2 = dx * dx + dy * dy + dz * dz;
          const real pulse = 0.01f * std::exp(-r2 / 6.0f);
          state.rho_d(i, j, k) += pulse;
          state.rho_theta_m(i, j, k) = state.rho_d(i, j, k) * 300.0f;
        }
      }
    }

    const double mass_before = state.total_dry_mass();

    dycore::DryStepperConfig cfg_step{};
    cfg_step.dt = 0.5f;
    cfg_step.fast_substeps = 6;
    dycore::LocalSplitExplicitFastMode fast_modes;
    for (int n = 0; n < 8; ++n) {
      fast_modes.apply(states, domain.layout, domain.metrics, cfg_step);
    }

    const auto summary = dycore::summarize_dry_states(states);
    TEST_NEAR(mass_before, state.total_dry_mass(), 1.0e-3);
    TEST_CHECK(summary.u_face.max > 1.0e-4);
    TEST_CHECK(summary.v_face.max > 1.0e-4);
    TEST_NEAR(state.rho_d(ic + 2, jc, kc), state.rho_d(ic - 2, jc, kc), 5.0e-4);
    TEST_NEAR(state.rho_d(ic, jc + 2, kc), state.rho_d(ic, jc - 2, kc), 5.0e-4);
    TEST_NEAR(summary.total_horizontal_momentum_x, 0.0, 1.0e-3);
    TEST_NEAR(summary.total_horizontal_momentum_y, 0.0, 1.0e-3);
    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
