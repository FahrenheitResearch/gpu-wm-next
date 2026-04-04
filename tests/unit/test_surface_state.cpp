#include <stdexcept>

#include "gwm/surface/surface_state.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm::surface;

  SurfaceState state(3, 2, 2, 4);
  TEST_CHECK(state.nx() == 3);
  TEST_CHECK(state.ny() == 2);
  TEST_CHECK(state.ntile() == 2);
  TEST_CHECK(state.nsoil() == 4);
  TEST_CHECK(state.cell_count() == 6);
  TEST_CHECK(state.tile_cell_count() == 12);
  TEST_CHECK(state.soil_storage_count() == 48);
  TEST_CHECK(state.tile_linear_index(1, 2, 1) == 11);
  TEST_CHECK(state.soil_linear_index(3, 1, 2, 1) == 47);

  state.skin_temperature(1, 2, 1) = 301.5f;
  state.canopy_water(0, 1, 0) = 0.2f;
  state.snow_water(1, 0, 1) = 0.7f;
  state.soil_temperature(3, 1, 2, 1) = 289.0f;
  state.soil_moisture(2, 0, 1, 0) = 0.33f;
  state.validate();

  const SurfaceState& const_state = state;
  TEST_NEAR(const_state.skin_temperature(1, 2, 1), 301.5f, 1.0e-6f);
  TEST_NEAR(const_state.canopy_water(0, 1, 0), 0.2f, 1.0e-6f);
  TEST_NEAR(const_state.snow_water(1, 0, 1), 0.7f, 1.0e-6f);
  TEST_NEAR(const_state.soil_temperature(3, 1, 2, 1), 289.0f, 1.0e-6f);
  TEST_NEAR(const_state.soil_moisture(2, 0, 1, 0), 0.33f, 1.0e-6f);

  SurfaceState single_tile(2, 2, 1, 4);
  TEST_CHECK(single_tile.ntile() == 1);
  TEST_CHECK(single_tile.tile_cell_count() == single_tile.cell_count());
  single_tile.validate();

  bool threw = false;
  try {
    SurfaceState bad(0, 2, 1, 4);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  TEST_CHECK(threw);

  threw = false;
  try {
    (void)state.skin_temperature(2, 0, 0);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  TEST_CHECK(threw);

  threw = false;
  try {
    (void)state.soil_temperature(4, 0, 0, 0);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  TEST_CHECK(threw);

  return 0;
}
