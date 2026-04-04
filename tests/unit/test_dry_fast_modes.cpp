#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/dry_diagnostics.hpp"

#include "test_assert.hpp"

namespace {

double abs_value(double value) { return value < 0.0 ? -value : value; }

double total_dry_mass(const std::vector<gwm::dycore::DryState>& states) {
  double total = 0.0;
  for (const auto& state : states) {
    total += state.total_dry_mass();
  }
  return total;
}

gwm::real max_abs_face_value(const gwm::dycore::ScalarRange& range) {
  const double abs_min = abs_value(range.min);
  const double abs_max = abs_value(range.max);
  return static_cast<gwm::real>(abs_min > abs_max ? abs_min : abs_max);
}

void seed_center_density_pulse(
    std::vector<gwm::dycore::DryState>& states,
    const std::vector<gwm::domain::SubdomainDescriptor>& layout,
    const gwm::domain::RectilinearDomainConfig& cfg, gwm::real rho_perturbation,
    gwm::real theta_ref) {
  for (std::size_t n = 0; n < states.size(); ++n) {
    auto& state = states[n];
    const auto& desc = layout[n];
    for (int k = 0; k < state.rho_d.nz(); ++k) {
      for (int j = 0; j < state.rho_d.ny(); ++j) {
        for (int i = 0; i < state.rho_d.nx(); ++i) {
          const int global_i = desc.i_begin + i;
          const int global_j = desc.j_begin + j;
          if (global_i == cfg.nx / 2 && global_j == cfg.ny / 2 &&
              k == cfg.nz / 2) {
            state.rho_d(i, j, k) += rho_perturbation;
            state.rho_theta_m(i, j, k) = state.rho_d(i, j, k) * theta_ref;
          }
        }
      }
    }
  }
}

}  // namespace

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
      seed_center_density_pulse(states, domain.layout, cfg, 0.02f, 300.0f);

      const double mass_before = total_dry_mass(states);

      dycore::DryStepperConfig step_cfg{};
      step_cfg.dt = 1.0f;
      step_cfg.fast_substeps = 4;
      dycore::LocalSplitExplicitFastMode fast_modes;
      for (int n = 0; n < 3; ++n) {
        fast_modes.apply(states, domain.layout, domain.metrics, step_cfg);
      }

      const double mass_after = total_dry_mass(states);
      TEST_NEAR(mass_before, mass_after, 1.0e-3);
    }

    {
      auto states_horizontal_on = dycore::make_constant_dry_state(
          domain.layout, 1.0f, 300.0f, 0.0f, 0.0f, 0.0f, "toggle_hmom_on");
      auto states_horizontal_off = dycore::make_constant_dry_state(
          domain.layout, 1.0f, 300.0f, 0.0f, 0.0f, 0.0f, "toggle_hmom_off");
      seed_center_density_pulse(states_horizontal_on, domain.layout, cfg, 0.02f,
                                300.0f);
      seed_center_density_pulse(states_horizontal_off, domain.layout, cfg,
                                0.02f, 300.0f);

      const double mass_before = total_dry_mass(states_horizontal_on);

      dycore::DryStepperConfig cfg_horizontal_on{};
      cfg_horizontal_on.dt = 1.0f;
      cfg_horizontal_on.fast_substeps = 1;
      cfg_horizontal_on.fast_update_horizontal_momentum = true;
      cfg_horizontal_on.fast_update_density = false;

      dycore::DryStepperConfig cfg_horizontal_off = cfg_horizontal_on;
      cfg_horizontal_off.fast_update_horizontal_momentum = false;

      dycore::LocalSplitExplicitFastMode fast_modes;
      fast_modes.apply(states_horizontal_on, domain.layout, domain.metrics,
                       cfg_horizontal_on);
      fast_modes.apply(states_horizontal_off, domain.layout, domain.metrics,
                       cfg_horizontal_off);

      const auto summary_horizontal_on =
          dycore::summarize_dry_states(states_horizontal_on);
      const auto summary_horizontal_off =
          dycore::summarize_dry_states(states_horizontal_off);
      const auto on_u_mag = max_abs_face_value(summary_horizontal_on.u_face);
      const auto on_v_mag = max_abs_face_value(summary_horizontal_on.v_face);
      const auto off_u_mag = max_abs_face_value(summary_horizontal_off.u_face);
      const auto off_v_mag = max_abs_face_value(summary_horizontal_off.v_face);

      TEST_NEAR(mass_before, total_dry_mass(states_horizontal_on), 1.0e-5);
      TEST_NEAR(mass_before, total_dry_mass(states_horizontal_off), 1.0e-5);
      TEST_CHECK(on_u_mag > 1.0e-5f);
      TEST_CHECK(on_v_mag > 1.0e-5f);
      TEST_CHECK(off_u_mag < 1.0e-8f);
      TEST_CHECK(off_v_mag < 1.0e-8f);
      TEST_NEAR(summary_horizontal_on.rho_d.max, 1.02, 1.0e-6);
      TEST_NEAR(summary_horizontal_off.rho_d.max, 1.02, 1.0e-6);
      TEST_NEAR(summary_horizontal_on.rho_d.min, 1.0, 1.0e-6);
      TEST_NEAR(summary_horizontal_off.rho_d.min, 1.0, 1.0e-6);
    }

    {
      auto states_density_on = dycore::make_constant_dry_state(
          domain.layout, 1.0f, 300.0f, 0.0f, 0.0f, 0.0f, "toggle_rho_on");
      auto states_density_off = dycore::make_constant_dry_state(
          domain.layout, 1.0f, 300.0f, 0.0f, 0.0f, 0.0f, "toggle_rho_off");
      seed_center_density_pulse(states_density_on, domain.layout, cfg, 0.02f,
                                300.0f);
      seed_center_density_pulse(states_density_off, domain.layout, cfg, 0.02f,
                                300.0f);

      const double mass_before = total_dry_mass(states_density_on);

      dycore::DryStepperConfig cfg_density_on{};
      cfg_density_on.dt = 1.0f;
      cfg_density_on.fast_substeps = 1;
      cfg_density_on.fast_update_horizontal_momentum = false;
      cfg_density_on.fast_update_density = true;

      dycore::DryStepperConfig cfg_density_off = cfg_density_on;
      cfg_density_off.fast_update_density = false;

      dycore::LocalSplitExplicitFastMode fast_modes;
      fast_modes.apply(states_density_on, domain.layout, domain.metrics,
                       cfg_density_on);
      fast_modes.apply(states_density_off, domain.layout, domain.metrics,
                       cfg_density_off);

      const auto summary_density_on =
          dycore::summarize_dry_states(states_density_on);
      const auto summary_density_off =
          dycore::summarize_dry_states(states_density_off);

      TEST_NEAR(mass_before, total_dry_mass(states_density_on), 1.0e-5);
      TEST_NEAR(mass_before, total_dry_mass(states_density_off), 1.0e-5);
      TEST_NEAR(summary_density_off.rho_d.max, 1.02, 1.0e-6);
      TEST_NEAR(summary_density_off.rho_d.min, 1.0, 1.0e-6);
      TEST_CHECK(abs_value(summary_density_on.rho_d.max - 1.02) > 1.0e-5);
    }

    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
