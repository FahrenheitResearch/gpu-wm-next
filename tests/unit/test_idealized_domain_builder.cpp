#include "gwm/domain/idealized_domain_builder.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm::domain;

  RectilinearDomainConfig config{};
  config.nx = 48;
  config.ny = 24;
  config.nz = 20;
  config.halo = 2;
  config.ranks_x = 2;
  config.ranks_y = 2;
  config.dx = 3000.0;
  config.dy = 3000.0;
  config.z_top = 18000.0;

  const auto domain = build_rectilinear_domain(config);
  TEST_CHECK(domain.metrics.nx == 48);
  TEST_CHECK(domain.metrics.nz == 20);
  TEST_NEAR(domain.metrics.dz_nominal, 900.0, 1.0e-9);
  TEST_CHECK(domain.layout.size() == 4);
  TEST_CHECK(domain.layout[0].halo == 2);
  TEST_CHECK(domain.layout[0].nx_local() == 24);
  TEST_CHECK(domain.layout[0].ny_local() == 12);
  return 0;
}
