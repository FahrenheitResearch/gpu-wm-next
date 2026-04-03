#include "gwm/domain/idealized_domain_builder.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm::domain;

  RectilinearDomainConfig config{};
  config.nx = 21;
  config.ny = 21;
  config.nz = 8;
  config.terrain_kind = TerrainProfileKind::CosineMountain;
  config.mountain_center_x = 0.5;
  config.mountain_center_y = 0.5;
  config.mountain_half_width_x = 0.25;
  config.mountain_half_width_y = 0.25;
  config.mountain_height = 1200.0;

  const auto domain = build_rectilinear_domain(config);
  TEST_NEAR(domain.terrain(10, 10), 1200.0f, 1.0e-3f);
  TEST_NEAR(domain.terrain(0, 0), 0.0f, 1.0e-6f);
  TEST_NEAR(domain.terrain(20, 20), 0.0f, 1.0e-6f);
  return 0;
}
