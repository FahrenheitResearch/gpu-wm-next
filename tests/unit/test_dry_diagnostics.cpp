#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_diagnostics.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;

  domain::RectilinearDomainConfig config{};
  config.nx = 4;
  config.ny = 3;
  config.nz = 2;
  config.halo = 1;

  const auto domain = domain::build_rectilinear_domain(config);
  auto states = dycore::make_constant_dry_state(domain.layout, 1.5f, 300.0f,
                                                2.0f, -1.0f, 0.25f,
                                                "diag_const");
  const auto summary = dycore::summarize_dry_states(states);

  TEST_NEAR(summary.total_dry_mass, 4.0 * 3.0 * 2.0 * 1.5, 1.0e-6);
  TEST_NEAR(summary.total_rho_theta_m, 4.0 * 3.0 * 2.0 * 1.5 * 300.0, 1.0e-6);
  TEST_NEAR(summary.total_horizontal_momentum_x,
            5.0 * 3.0 * 2.0 * 1.5 * 2.0, 1.0e-6);
  TEST_NEAR(summary.total_horizontal_momentum_y,
            4.0 * 4.0 * 2.0 * 1.5 * -1.0, 1.0e-6);
  TEST_NEAR(summary.total_vertical_momentum,
            4.0 * 3.0 * 3.0 * 1.5 * 0.25, 1.0e-6);
  TEST_NEAR(summary.rho_d.min, 1.5, 1.0e-6);
  TEST_NEAR(summary.rho_d.max, 1.5, 1.0e-6);
  TEST_NEAR(summary.theta_m.min, 300.0, 1.0e-6);
  TEST_NEAR(summary.theta_m.max, 300.0, 1.0e-6);
  TEST_NEAR(summary.u_face.min, 1.5 * 2.0, 1.0e-6);
  TEST_NEAR(summary.u_face.max, 1.5 * 2.0, 1.0e-6);
  TEST_NEAR(summary.v_face.min, 1.5 * -1.0, 1.0e-6);
  TEST_NEAR(summary.v_face.max, 1.5 * -1.0, 1.0e-6);
  TEST_NEAR(summary.w_face.min, 1.5 * 0.25, 1.0e-6);
  TEST_NEAR(summary.w_face.max, 1.5 * 0.25, 1.0e-6);

  return 0;
}
