#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_pressure_gradient.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;

  domain::RectilinearDomainConfig cfg{};
  cfg.nx = 6;
  cfg.ny = 5;
  cfg.nz = 2;
  cfg.halo = 1;
  cfg.ranks_x = 1;
  cfg.ranks_y = 1;
  cfg.periodic_x = false;
  cfg.periodic_y = false;

  const auto domain = domain::build_rectilinear_domain(cfg);

  std::vector<state::Field3D<real>> pressure_fields;
  pressure_fields.emplace_back(cfg.nx, cfg.ny, cfg.nz, cfg.halo, "pressure");
  pressure_fields[0].fill(100000.0f);

  std::vector<dycore::DrySlowTendencies> tendencies;
  tendencies.emplace_back(cfg.nx, cfg.ny, cfg.nz, cfg.halo, "tendency");
  tendencies[0].fill_zero();

  dycore::add_horizontal_pressure_gradient_tendencies(
      pressure_fields, domain.layout, domain.metrics, tendencies);

  TEST_NEAR(tendencies[0].mom_u.storage()(1, 2, 0), 0.0f, 1.0e-6f);
  TEST_NEAR(tendencies[0].mom_v.storage()(2, 1, 0), 0.0f, 1.0e-6f);

  for (int k = 0; k < cfg.nz; ++k) {
    for (int j = -cfg.halo; j < cfg.ny + cfg.halo; ++j) {
      for (int i = -cfg.halo; i < cfg.nx + cfg.halo; ++i) {
        pressure_fields[0](i, j, k) =
            static_cast<real>(100000.0 + 10.0 * i + 20.0 * j);
      }
    }
  }
  tendencies[0].fill_zero();
  dycore::add_horizontal_pressure_gradient_tendencies(
      pressure_fields, domain.layout, domain.metrics, tendencies);

  TEST_NEAR(tendencies[0].mom_u.storage()(0, 2, 0), 0.0f, 1.0e-6f);
  TEST_NEAR(tendencies[0].mom_u.storage()(1, 2, 0),
            static_cast<real>(-10.0 / domain.metrics.dx), 1.0e-6f);
  TEST_NEAR(tendencies[0].mom_u.storage()(cfg.nx, 2, 0), 0.0f, 1.0e-6f);

  TEST_NEAR(tendencies[0].mom_v.storage()(2, 0, 0), 0.0f, 1.0e-6f);
  TEST_NEAR(tendencies[0].mom_v.storage()(2, 1, 0),
            static_cast<real>(-20.0 / domain.metrics.dy), 1.0e-6f);
  TEST_NEAR(tendencies[0].mom_v.storage()(2, cfg.ny, 0), 0.0f, 1.0e-6f);

  return 0;
}
