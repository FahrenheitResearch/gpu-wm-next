#include <stdexcept>

#include "gwm/ingest/canonical_ir.hpp"
#include "gwm/surface/surface_runtime_init.hpp"
#include "gwm/surface/surface_static_properties.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm::surface;

  SurfaceStaticProperties props(2, 2, 2);
  TEST_CHECK(props.nx() == 2);
  TEST_CHECK(props.ny() == 2);
  TEST_CHECK(props.ntile() == 2);
  TEST_CHECK(!props.is_single_tile());
  TEST_CHECK(props.cell_count() == 4);
  TEST_CHECK(props.tile_cell_count() == 8);
  TEST_CHECK(props.tile_linear_index(1, 1, 1) == 7);

  for (int j = 0; j < props.ny(); ++j) {
    for (int i = 0; i < props.nx(); ++i) {
      props.set_tile_fraction(0, i, j, 0.6f);
      props.set_tile_fraction(1, i, j, 0.4f);
      props.set_z0m(0, i, j, 0.08f);
      props.set_z0m(1, i, j, 0.2f);
      props.set_z0h(0, i, j, 0.02f);
      props.set_z0h(1, i, j, 0.05f);
      props.set_kind(0, i, j, SurfaceKind::Land);
      props.set_kind(1, i, j, SurfaceKind::Water);
    }
  }
  props.validate();

  const SurfaceStaticProperties& const_props = props;
  TEST_NEAR(const_props.tile_fraction(0, 1, 1), 0.6f, 1.0e-6f);
  TEST_NEAR(const_props.tile_fraction(1, 1, 1), 0.4f, 1.0e-6f);
  TEST_NEAR(const_props.z0m(1, 0, 0), 0.2f, 1.0e-6f);
  TEST_NEAR(const_props.z0h(0, 0, 0), 0.02f, 1.0e-6f);
  TEST_NEAR(const_props.tile_fraction_sum(0, 0), 1.0f, 1.0e-6f);
  TEST_CHECK(const_props.kind(1, 0, 0) == SurfaceKind::Water);

  SurfaceStaticProperties single_tile(1, 1, 1);
  TEST_CHECK(single_tile.is_single_tile());
  single_tile.validate();
  TEST_NEAR(single_tile.tile_fraction_sum(0, 0), 1.0f, 1.0e-6f);

  bool bad_fraction_threw = false;
  try {
    SurfaceStaticProperties bad_props(1, 1, 2);
    bad_props.set_tile_fraction(0, 0, 0, 0.7f);
    bad_props.set_tile_fraction(1, 0, 0, 0.1f);
    bad_props.validate();
  } catch (const std::runtime_error&) {
    bad_fraction_threw = true;
  }
  TEST_CHECK(bad_fraction_threw);

  bool bad_roughness_threw = false;
  try {
    SurfaceStaticProperties bad_props(1, 1, 1);
    bad_props.set_z0m(0, 0, 0, 0.0f);
  } catch (const std::runtime_error&) {
    bad_roughness_threw = true;
  }
  TEST_CHECK(bad_roughness_threw);

  bool bad_index_threw = false;
  try {
    (void)props.z0h(2, 0, 0);
  } catch (const std::runtime_error&) {
    bad_index_threw = true;
  }
  TEST_CHECK(bad_index_threw);

  gwm::ingest::FieldBundle2D static_surface_fields{};
  static_surface_fields.values["terrain_height"] = {10.0f, 20.0f, 30.0f, 40.0f};
  static_surface_fields.values["land_mask"] = {1.0f, 0.0f, 1.0f, 0.0f};
  static_surface_fields.values["land_use_index"] = {1.0f, 2.0f, 3.0f, 4.0f};

  const auto initialized =
      make_surface_static_properties_from_canonical_fields(static_surface_fields,
                                                           2, 2, 1);
  TEST_CHECK(initialized.is_single_tile());
  TEST_NEAR(initialized.tile_fraction_sum(0, 0), 1.0f, 1.0e-6f);
  TEST_CHECK(initialized.kind(0, 0, 0) == SurfaceKind::Land);
  TEST_CHECK(initialized.kind(0, 1, 0) == SurfaceKind::Water);
  TEST_NEAR(initialized.z0m(0, 0, 0), 0.10f, 1.0e-6f);
  TEST_NEAR(initialized.z0m(0, 1, 0), 0.0002f, 1.0e-6f);
  TEST_NEAR(initialized.z0h(0, 0, 0), 0.02f, 1.0e-6f);
  TEST_NEAR(initialized.z0h(0, 1, 0), 0.0002f, 1.0e-6f);
  initialized.validate();

  bool missing_land_mask_threw = false;
  try {
    gwm::ingest::FieldBundle2D missing_land_mask{};
    missing_land_mask.values["terrain_height"] = {10.0f, 20.0f, 30.0f, 40.0f};
    missing_land_mask.values["land_use_index"] = {1.0f, 2.0f, 3.0f, 4.0f};
    (void)make_surface_static_properties_from_canonical_fields(
        missing_land_mask, 2, 2, 1);
  } catch (const std::runtime_error&) {
    missing_land_mask_threw = true;
  }
  TEST_CHECK(missing_land_mask_threw);

  bool multi_tile_threw = false;
  try {
    (void)make_surface_static_properties_from_canonical_fields(
        static_surface_fields, 2, 2, 2);
  } catch (const std::runtime_error&) {
    multi_tile_threw = true;
  }
  TEST_CHECK(multi_tile_threw);

  return 0;
}
