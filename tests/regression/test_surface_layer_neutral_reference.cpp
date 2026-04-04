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

  TEST_NEAR(result.tile_outputs[0].t2, 301.46582f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[0].q2, 0.010266f, 1.0e-5f);
  TEST_NEAR(result.tile_outputs[0].u10, 8.63242f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[0].v10, 2.87747f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[0].cm, 0.009800f, 1.0e-5f);
  TEST_NEAR(result.tile_outputs[0].ustar, 1.24097f, 1.0e-4f);

  TEST_NEAR(result.tile_outputs[1].t2, 302.88330f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[1].q2, 0.012117f, 1.0e-5f);
  TEST_NEAR(result.tile_outputs[1].u10, 7.41112f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[1].v10, 2.47037f, 1.0e-4f);
  TEST_NEAR(result.tile_outputs[1].cm, 0.016686f, 1.0e-5f);
  TEST_NEAR(result.tile_outputs[1].ustar, 1.62039f, 1.0e-4f);

  TEST_NEAR(result.cell_mean.t2, 301.89105f, 1.0e-4f);
  TEST_NEAR(result.cell_mean.q2, 0.010821f, 1.0e-5f);
  TEST_NEAR(result.cell_mean.u10, 8.26603f, 1.0e-4f);
  TEST_NEAR(result.cell_mean.v10, 2.75534f, 1.0e-4f);
  TEST_NEAR(result.cell_mean.cm, 0.011866f, 1.0e-5f);
  TEST_NEAR(result.cell_mean.ch, 0.008386f, 1.0e-5f);
  TEST_NEAR(result.cell_mean.cq, 0.008386f, 1.0e-5f);
  TEST_NEAR(result.cell_mean.ustar, 1.35480f, 1.0e-4f);

  return 0;
}
