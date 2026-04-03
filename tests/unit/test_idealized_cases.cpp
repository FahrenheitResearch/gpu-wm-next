#include <algorithm>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/idealized_cases.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;

  domain::RectilinearDomainConfig config{};
  config.nx = 32;
  config.ny = 32;
  config.nz = 16;
  config.ranks_x = 2;
  config.ranks_y = 2;
  config.dx = 1000.0;
  config.dy = 1000.0;
  config.z_top = 16000.0;

  const auto domain = domain::build_rectilinear_domain(config);

  dycore::ThermoBubbleConfig bubble_cfg{};
  bubble_cfg.rho_background = 1.0f;
  bubble_cfg.theta_background = 300.0f;
  bubble_cfg.u_background = 10.0f;
  bubble_cfg.theta_perturbation = 3.0f;
  bubble_cfg.center_x_fraction = 0.5f;
  bubble_cfg.center_y_fraction = 0.5f;
  bubble_cfg.center_z_fraction = 0.4f;
  bubble_cfg.radius_x_fraction = 0.12f;
  bubble_cfg.radius_y_fraction = 0.12f;
  bubble_cfg.radius_z_fraction = 0.15f;

  auto warm = dycore::make_warm_bubble_state(domain.layout, domain.metrics,
                                             bubble_cfg, "warm");
  auto cold = dycore::make_density_current_state(domain.layout, domain.metrics,
                                                 bubble_cfg, "cold");

  const auto warm_theta =
      comm::VirtualRankLayout::gather_scalar(
          [&warm]() {
            std::vector<state::Field3D<real>> fields;
            fields.reserve(warm.size());
            for (const auto& state : warm) {
              fields.push_back(state.rho_theta_m.clone_empty_like("_gather"));
              fields.back().copy_all_from(state.rho_theta_m);
            }
            return fields;
          }(),
          domain.layout, "warm_theta");

  const auto cold_theta =
      comm::VirtualRankLayout::gather_scalar(
          [&cold]() {
            std::vector<state::Field3D<real>> fields;
            fields.reserve(cold.size());
            for (const auto& state : cold) {
              fields.push_back(state.rho_theta_m.clone_empty_like("_gather"));
              fields.back().copy_all_from(state.rho_theta_m);
            }
            return fields;
          }(),
          domain.layout, "cold_theta");

  real warm_max = -1.0e9f;
  real cold_min = 1.0e9f;
  for (int k = 0; k < warm_theta.nz(); ++k) {
    for (int j = 0; j < warm_theta.ny(); ++j) {
      for (int i = 0; i < warm_theta.nx(); ++i) {
        warm_max = std::max(warm_max, warm_theta(i, j, k));
        cold_min = std::min(cold_min, cold_theta(i, j, k));
      }
    }
  }

  TEST_CHECK(warm_max > bubble_cfg.rho_background * bubble_cfg.theta_background);
  TEST_CHECK(cold_min < bubble_cfg.rho_background * bubble_cfg.theta_background);
  TEST_CHECK(cold_min > 0.0f);

  return 0;
}
