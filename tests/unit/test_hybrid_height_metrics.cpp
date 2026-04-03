#include "gwm/domain/idealized_domain_builder.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm::domain;

  {
    RectilinearDomainConfig flat{};
    flat.nx = 12;
    flat.ny = 10;
    flat.nz = 8;
    flat.z_top = 16000.0;

    const auto domain = build_rectilinear_domain(flat);
    TEST_NEAR(domain.metrics.z_interface(0, 0, 0), 0.0f, 1.0e-6f);
    TEST_NEAR(domain.metrics.z_interface(0, 0, flat.nz),
              static_cast<gwm::real>(flat.z_top), 1.0e-4f);
    TEST_NEAR(domain.metrics.inv_dz_cell(0, 0, 0),
              static_cast<gwm::real>(1.0 / domain.metrics.dz_nominal), 1.0e-7f);
    for (int k = 0; k < flat.nz; ++k) {
      TEST_NEAR(domain.metrics.z_center(0, 0, k),
                static_cast<gwm::real>(domain.metrics.flat_z_center(k)), 1.0e-4f);
      TEST_CHECK(domain.metrics.z_interface(0, 0, k + 1) >
                 domain.metrics.z_interface(0, 0, k));
    }
  }

  {
    RectilinearDomainConfig mountain{};
    mountain.nx = 21;
    mountain.ny = 21;
    mountain.nz = 12;
    mountain.z_top = 18000.0;
    mountain.terrain_kind = TerrainProfileKind::CosineMountain;
    mountain.mountain_center_x = 0.5;
    mountain.mountain_center_y = 0.5;
    mountain.mountain_half_width_x = 0.25;
    mountain.mountain_half_width_y = 0.25;
    mountain.mountain_height = 1200.0;
    mountain.terrain_taper_eta = 0.25;

    const auto domain = build_rectilinear_domain(mountain);
    const int ic = mountain.nx / 2;
    const int jc = mountain.ny / 2;

    TEST_NEAR(domain.metrics.z_interface(ic, jc, 0), domain.terrain(ic, jc), 1.0e-3f);
    TEST_NEAR(domain.metrics.z_interface(ic, jc, mountain.nz),
              static_cast<gwm::real>(mountain.z_top), 1.0e-4f);

    for (int k = 0; k < mountain.nz; ++k) {
      TEST_CHECK(domain.metrics.z_interface(ic, jc, k + 1) >
                 domain.metrics.z_interface(ic, jc, k));
      TEST_CHECK(domain.metrics.z_center(ic, jc, k) >
                 domain.metrics.z_interface(ic, jc, k));
      TEST_CHECK(domain.metrics.z_center(ic, jc, k) <
                 domain.metrics.z_interface(ic, jc, k + 1));
      TEST_CHECK(domain.metrics.inv_dz_cell(ic, jc, k) > 0.0f);
    }

    const auto low_terrain_offset =
        domain.metrics.z_center(ic, jc, 0) -
        static_cast<gwm::real>(domain.metrics.flat_z_center(0));
    const auto upper_terrain_offset =
        domain.metrics.z_center(ic, jc, mountain.nz - 1) -
        static_cast<gwm::real>(domain.metrics.flat_z_center(mountain.nz - 1));
    TEST_CHECK(low_terrain_offset > upper_terrain_offset);
    TEST_CHECK(domain.metrics.terrain_weight_center(mountain.nz - 1) < 5.0e-2);
  }

  return 0;
}
