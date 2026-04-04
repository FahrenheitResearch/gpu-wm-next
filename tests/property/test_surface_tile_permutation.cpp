#include "gwm/surface/surface_layer_exchange.hpp"

#include "test_assert.hpp"

namespace {

void configure_two_tile_case(gwm::surface::SurfaceState& state,
                             gwm::surface::SurfaceStaticProperties& props) {
  state.skin_temperature(0, 0, 0) = 301.0f;
  state.skin_temperature(1, 0, 0) = 297.0f;

  props.set_tile_fraction(0, 0, 0, 0.65f);
  props.set_tile_fraction(1, 0, 0, 0.35f);
  props.set_z0m(0, 0, 0, 0.08f);
  props.set_z0m(1, 0, 0, 0.30f);
  props.set_z0h(0, 0, 0, 0.02f);
  props.set_z0h(1, 0, 0, 0.06f);
}

void check_diag_near(const gwm::surface::SurfaceLayerDiagnostics& a,
                     const gwm::surface::SurfaceLayerDiagnostics& b,
                     gwm::real tol) {
  TEST_NEAR(a.t2, b.t2, tol);
  TEST_NEAR(a.q2, b.q2, tol);
  TEST_NEAR(a.rh2, b.rh2, tol);
  TEST_NEAR(a.u10, b.u10, tol);
  TEST_NEAR(a.v10, b.v10, tol);
  TEST_NEAR(a.wind_speed_ref, b.wind_speed_ref, tol);
  TEST_NEAR(a.wind_speed_target, b.wind_speed_target, tol);
  TEST_NEAR(a.cm, b.cm, tol);
  TEST_NEAR(a.ch, b.ch, tol);
  TEST_NEAR(a.cq, b.cq, tol);
  TEST_NEAR(a.ustar, b.ustar, tol);
}

}  // namespace

int main() {
  using namespace gwm;

  surface::SurfaceLayerExchangeForcing forcing{};
  forcing.z_ref = 35.0f;
  forcing.z_target_temp = 2.0f;
  forcing.z_target_wind = 10.0f;
  forcing.u_ref = 14.0f;
  forcing.v_ref = -3.0f;
  forcing.theta_ref = 303.0f;
  forcing.q_ref = 0.009f;
  forcing.psfc = 94000.0f;
  forcing.tile_q_surface = {0.0125f, 0.0080f};

  surface::SurfaceState state_a(1, 1, 2, 4);
  surface::SurfaceStaticProperties props_a(1, 1, 2);
  configure_two_tile_case(state_a, props_a);

  surface::SurfaceState state_b(1, 1, 2, 4);
  surface::SurfaceStaticProperties props_b(1, 1, 2);
  state_b.skin_temperature(0, 0, 0) = state_a.skin_temperature(1, 0, 0);
  state_b.skin_temperature(1, 0, 0) = state_a.skin_temperature(0, 0, 0);
  props_b.set_tile_fraction(0, 0, 0, props_a.tile_fraction(1, 0, 0));
  props_b.set_tile_fraction(1, 0, 0, props_a.tile_fraction(0, 0, 0));
  props_b.set_z0m(0, 0, 0, props_a.z0m(1, 0, 0));
  props_b.set_z0m(1, 0, 0, props_a.z0m(0, 0, 0));
  props_b.set_z0h(0, 0, 0, props_a.z0h(1, 0, 0));
  props_b.set_z0h(1, 0, 0, props_a.z0h(0, 0, 0));

  auto forcing_b = forcing;
  forcing_b.tile_q_surface = {forcing.tile_q_surface[1], forcing.tile_q_surface[0]};

  const auto result_a =
      surface::evaluate_neutral_surface_exchange(state_a, props_a, 0, 0, forcing);
  const auto result_b = surface::evaluate_neutral_surface_exchange(
      state_b, props_b, 0, 0, forcing_b);

  check_diag_near(result_a.cell_mean, result_b.cell_mean, 1.0e-6f);
  check_diag_near(result_a.tile_outputs[0], result_b.tile_outputs[1], 1.0e-6f);
  check_diag_near(result_a.tile_outputs[1], result_b.tile_outputs[0], 1.0e-6f);
  return 0;
}
