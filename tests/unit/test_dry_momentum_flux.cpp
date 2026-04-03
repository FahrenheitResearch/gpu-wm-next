#include <cmath>

#include "gwm/comm/halo_exchange.hpp"
#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/dry_momentum_flux.hpp"
#include "gwm/dycore/idealized_cases.hpp"

#include "test_assert.hpp"

namespace {

double face_owned_sum(const gwm::state::FaceField<gwm::real>& field) {
  double sum = 0.0;
  const auto& storage = field.storage();
  for (int k = 0; k < storage.nz(); ++k) {
    for (int j = 0; j < storage.ny(); ++j) {
      for (int i = 0; i < storage.nx(); ++i) {
        sum += static_cast<double>(storage(i, j, k));
      }
    }
  }
  return sum;
}

}  // namespace

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig cfg{};
    cfg.nx = 20;
    cfg.ny = 18;
    cfg.nz = 8;
    cfg.halo = 1;
    cfg.ranks_x = 2;
    cfg.ranks_y = 2;

    const auto domain = domain::build_rectilinear_domain(cfg);

    {
      auto states = dycore::make_constant_dry_state(
          domain.layout, 1.05f, 300.0f, 12.0f, -4.0f, 0.5f, "const_flux");

      std::vector<state::Field3D<real>> rho_fields;
      std::vector<state::FaceField<real>> mom_u_fields;
      std::vector<state::FaceField<real>> mom_v_fields;
      std::vector<state::FaceField<real>> mom_w_fields;
      std::vector<dycore::DrySlowTendencies> tendencies;
      for (const auto& state : states) {
        rho_fields.push_back(state.rho_d.clone_empty_like("_rho"));
        rho_fields.back().copy_all_from(state.rho_d);
        mom_u_fields.emplace_back(state.rho_d.nx(), state.rho_d.ny(),
                                  state.rho_d.nz(), state.rho_d.halo(),
                                  state::FaceOrientation::X, "u");
        mom_u_fields.back().storage().copy_all_from(state.mom_u.storage());
        mom_v_fields.emplace_back(state.rho_d.nx(), state.rho_d.ny(),
                                  state.rho_d.nz(), state.rho_d.halo(),
                                  state::FaceOrientation::Y, "v");
        mom_v_fields.back().storage().copy_all_from(state.mom_v.storage());
        mom_w_fields.emplace_back(state.rho_d.nx(), state.rho_d.ny(),
                                  state.rho_d.nz(), state.rho_d.halo(),
                                  state::FaceOrientation::Z, "w");
        mom_w_fields.back().storage().copy_all_from(state.mom_w.storage());
        tendencies.emplace_back(state.rho_d.nx(), state.rho_d.ny(),
                                state.rho_d.nz(), state.rho_d.halo(), "t");
        tendencies.back().fill_zero();
      }

      comm::HaloExchange::exchange_scalar(rho_fields, domain.layout);
      comm::HaloExchange::exchange_face(mom_u_fields, domain.layout);
      comm::HaloExchange::exchange_face(mom_v_fields, domain.layout);
      comm::HaloExchange::exchange_face(mom_w_fields, domain.layout);

      dycore::add_dry_momentum_flux_tendencies(
          rho_fields, mom_u_fields, mom_v_fields, mom_w_fields, domain.layout,
          domain.metrics, tendencies);

      for (const auto& tendency : tendencies) {
        TEST_NEAR(face_owned_sum(tendency.mom_u), 0.0, 1.0e-5);
        TEST_NEAR(face_owned_sum(tendency.mom_v), 0.0, 1.0e-5);
        TEST_NEAR(face_owned_sum(tendency.mom_w), 0.0, 1.0e-5);
      }
    }

    {
      dycore::ThermoBubbleConfig cfg_bubble{};
      cfg_bubble.rho_background = 1.0f;
      cfg_bubble.theta_background = 300.0f;
      cfg_bubble.theta_perturbation = 3.0f;
      cfg_bubble.center_x_fraction = 0.5f;
      cfg_bubble.center_y_fraction = 0.5f;
      cfg_bubble.center_z_fraction = 0.35f;
      cfg_bubble.radius_x_fraction = 0.15f;
      cfg_bubble.radius_y_fraction = 0.15f;
      cfg_bubble.radius_z_fraction = 0.2f;

      auto states = dycore::make_density_current_state(
          domain.layout, domain.metrics, cfg_bubble, "pulse_flux");
      dycore::DryStepperConfig step_cfg{};
      step_cfg.dt = 2.0f;
      step_cfg.fast_substeps = 3;
      dycore::NullBoundaryUpdater boundary;
      dycore::LocalSplitExplicitFastMode fast_modes;
      dycore::advance_dry_state_ssprk3(states, domain.layout, domain.metrics,
                                       step_cfg, boundary, fast_modes);

      bool any_nonzero_u = false;
      for (const auto& state : states) {
        const auto& storage = state.mom_u.storage();
        for (int k = 0; k < storage.nz(); ++k) {
          for (int j = 0; j < storage.ny(); ++j) {
            for (int i = 0; i < storage.nx(); ++i) {
              const real value = storage(i, j, k);
              TEST_CHECK(std::isfinite(static_cast<double>(value)));
              if (std::fabs(static_cast<double>(value)) > 1.0e-4) {
                any_nonzero_u = true;
              }
            }
          }
        }
      }
      TEST_CHECK(any_nonzero_u);
    }

    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
