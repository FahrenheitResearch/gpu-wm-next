#include "gwm/surface/surface_layer_exchange.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;

  surface::SurfaceState state(1, 1, 2, 4);
  state.skin_temperature(0, 0, 0) = 300.0f;
  state.skin_temperature(1, 0, 0) = 302.0f;

  surface::SurfaceStaticProperties props(1, 1, 2);
  props.set_tile_fraction(0, 0, 0, 0.7f);
  props.set_tile_fraction(1, 0, 0, 0.3f);
  props.set_z0m(0, 0, 0, 0.1f);
  props.set_z0m(1, 0, 0, 0.3f);
  props.set_z0h(0, 0, 0, 0.02f);
  props.set_z0h(1, 0, 0, 0.06f);

  surface::SurfaceLayerExchangeForcing forcing{};
  forcing.z_ref = 30.0f;
  forcing.z_target_temp = 2.0f;
  forcing.z_target_wind = 10.0f;
  forcing.u_ref = 12.0f;
  forcing.v_ref = 4.0f;
  forcing.theta_ref = 304.0f;
  forcing.q_ref = 0.009f;
  forcing.psfc = 96000.0f;
  forcing.tile_q_surface = {0.011f, 0.013f};

  const auto result =
      surface::evaluate_neutral_surface_exchange(state, props, 0, 0, forcing);

  TEST_NEAR(result.tile_outputs[0].t2, 302.51883f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[0].q2, 0.009741f, 1.0e-5f);
  TEST_NEAR(result.tile_outputs[0].u10, 9.68867f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[0].v10, 3.22956f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[0].cm, 0.004918f, 1.0e-5f);
  TEST_NEAR(result.tile_outputs[0].ustar, 0.88707f, 1.0e-4f);

  TEST_NEAR(result.tile_outputs[1].t2, 303.12848f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[1].q2, 0.010743f, 1.0e-5f);
  TEST_NEAR(result.tile_outputs[1].u10, 9.13727f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[1].v10, 3.04576f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[1].cm, 0.007544f, 1.0e-5f);
  TEST_NEAR(result.tile_outputs[1].ustar, 1.09869f, 1.0e-4f);

  TEST_NEAR(result.cell_mean.t2, 302.70172f, 1.0e-4f);
  TEST_NEAR(result.cell_mean.q2, 0.010041f, 1.0e-5f);
  TEST_NEAR(result.cell_mean.u10, 9.52325f, 1.0e-4f);
  TEST_NEAR(result.cell_mean.v10, 3.17442f, 1.0e-4f);
  TEST_NEAR(result.cell_mean.cm, 0.005706f, 1.0e-5f);
  TEST_NEAR(result.cell_mean.ch, 0.004362f, 1.0e-5f);
  TEST_NEAR(result.cell_mean.cq, 0.004362f, 1.0e-5f);
  TEST_NEAR(result.cell_mean.ustar, 0.95055f, 1.0e-4f);

  return 0;
}
