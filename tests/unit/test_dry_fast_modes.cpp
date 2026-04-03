#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/dry_diagnostics.hpp"

#include "test_assert.hpp"

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig cfg{};
    cfg.nx = 20;
    cfg.ny = 20;
    cfg.nz = 8;
    cfg.halo = 1;
    cfg.ranks_x = 2;
    cfg.ranks_y = 2;

    const auto domain = domain::build_rectilinear_domain(cfg);

    {
      auto states = dycore::make_constant_dry_state(
          domain.layout, 1.1f, 300.0f, 12.0f, -6.0f, 0.0f, "const_fast");
      const auto before = dycore::summarize_dry_states(states);

      dycore::DryStepperConfig step_cfg{};
      step_cfg.dt = 2.0f;
      step_cfg.fast_substeps = 4;
      dycore::LocalSplitExplicitFastMode fast_modes;
      fast_modes.apply(states, domain.layout, domain.metrics, step_cfg);

      const auto after = dycore::summarize_dry_states(states);
      TEST_NEAR(before.total_dry_mass, after.total_dry_mass, 1.0e-4);
      TEST_NEAR(after.total_horizontal_momentum_x,
                before.total_horizontal_momentum_x, 1.0e-4);
      TEST_NEAR(after.total_horizontal_momentum_y,
                before.total_horizontal_momentum_y, 1.0e-4);
      TEST_NEAR(after.u_face.min, before.u_face.min, 1.0e-5);
      TEST_NEAR(after.u_face.max, before.u_face.max, 1.0e-5);
      TEST_NEAR(after.v_face.min, before.v_face.min, 1.0e-5);
      TEST_NEAR(after.v_face.max, before.v_face.max, 1.0e-5);
    }

    {
      auto states = dycore::make_constant_dry_state(
          domain.layout, 1.0f, 300.0f, 0.0f, 0.0f, 0.0f, "pulse_fast");
      for (std::size_t n = 0; n < states.size(); ++n) {
        auto& state = states[n];
        const auto& desc = domain.layout[n];
        for (int k = 0; k < state.rho_d.nz(); ++k) {
          for (int j = 0; j < state.rho_d.ny(); ++j) {
            for (int i = 0; i < state.rho_d.nx(); ++i) {
              const int global_i = desc.i_begin + i;
              const int global_j = desc.j_begin + j;
              if (global_i == cfg.nx / 2 && global_j == cfg.ny / 2 &&
                  k == cfg.nz / 2) {
                state.rho_d(i, j, k) = 1.02f;
                state.rho_theta_m(i, j, k) = 1.02f * 300.0f;
              }
            }
          }
        }
      }

      double mass_before = 0.0;
      for (const auto& state : states) {
        mass_before += state.total_dry_mass();
      }

      dycore::DryStepperConfig step_cfg{};
      step_cfg.dt = 1.0f;
      step_cfg.fast_substeps = 4;
      dycore::LocalSplitExplicitFastMode fast_modes;
      for (int n = 0; n < 3; ++n) {
        fast_modes.apply(states, domain.layout, domain.metrics, step_cfg);
      }

      double mass_after = 0.0;
      for (const auto& state : states) {
        mass_after += state.total_dry_mass();
      }
      TEST_NEAR(mass_before, mass_after, 1.0e-3);
    }

    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
