#include <algorithm>
#include <limits>

#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/idealized_cases.hpp"

#include "test_assert.hpp"

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig domain_cfg{};
    domain_cfg.nx = 32;
    domain_cfg.ny = 32;
    domain_cfg.nz = 16;
    domain_cfg.halo = 1;
    domain_cfg.ranks_x = 2;
    domain_cfg.ranks_y = 2;
    domain_cfg.dx = 1000.0;
    domain_cfg.dy = 1000.0;
    domain_cfg.z_top = 16000.0;

    const auto domain = domain::build_rectilinear_domain(domain_cfg);

    dycore::ThermoBubbleConfig bubble_cfg{};
    bubble_cfg.rho_background = 1.0f;
    bubble_cfg.theta_background = 300.0f;
    bubble_cfg.theta_perturbation = 3.0f;
    bubble_cfg.center_x_fraction = 0.5f;
    bubble_cfg.center_y_fraction = 0.5f;
    bubble_cfg.center_z_fraction = 0.35f;
    bubble_cfg.radius_x_fraction = 0.12f;
    bubble_cfg.radius_y_fraction = 0.12f;
    bubble_cfg.radius_z_fraction = 0.16f;

    auto w_face_extrema = [](const std::vector<dycore::DryState>& states) {
      real max_value = -std::numeric_limits<real>::max();
      real min_value = std::numeric_limits<real>::max();
      for (const auto& state : states) {
        const auto& w = state.mom_w.storage();
        for (int k = 0; k < w.nz(); ++k) {
          for (int j = 0; j < w.ny(); ++j) {
            for (int i = 0; i < w.nx(); ++i) {
              const real value = w(i, j, k);
              max_value = std::max(max_value, value);
              min_value = std::min(min_value, value);
            }
          }
        }
      }
      return std::pair<real, real>{max_value, min_value};
    };

    dycore::DryStepperConfig step_cfg{};
    step_cfg.dt = 2.0f;
    step_cfg.fast_substeps = 3;
    dycore::NullBoundaryUpdater boundary;
    dycore::LocalSplitExplicitFastMode fast_modes;

    auto warm = dycore::make_warm_bubble_state(domain.layout, domain.metrics,
                                               bubble_cfg, "warm");
    auto [warm_initial_max, warm_initial_min] = w_face_extrema(warm);
    TEST_NEAR(warm_initial_max, 0.0f, 1.0e-6f);
    TEST_NEAR(warm_initial_min, 0.0f, 1.0e-6f);
    dycore::advance_dry_state_ssprk3(warm, domain.layout, domain.metrics,
                                     step_cfg, boundary, fast_modes);
    auto [warm_max, warm_min] = w_face_extrema(warm);
    TEST_CHECK(warm_max > 1.0e-3f);

    auto cold = dycore::make_density_current_state(domain.layout, domain.metrics,
                                                   bubble_cfg, "cold");
    auto [cold_initial_max, cold_initial_min] = w_face_extrema(cold);
    TEST_NEAR(cold_initial_max, 0.0f, 1.0e-6f);
    TEST_NEAR(cold_initial_min, 0.0f, 1.0e-6f);
    dycore::advance_dry_state_ssprk3(cold, domain.layout, domain.metrics,
                                     step_cfg, boundary, fast_modes);
    auto [cold_max, cold_min] = w_face_extrema(cold);
    TEST_CHECK(cold_min < -1.0e-3f);

    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
