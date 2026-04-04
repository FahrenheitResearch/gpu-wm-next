#include <algorithm>
#include <cmath>

#include "gwm/surface/surface_layer_closure.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;

  surface::SurfaceLayerInputs in{};
  in.z_ref = 40.0f;
  in.z_target_temp = 2.0f;
  in.z_target_wind = 10.0f;
  in.z0m = 0.1f;
  in.z0h = 0.02f;
  in.u_ref = 18.0f;
  in.v_ref = -6.0f;
  in.theta_ref = 304.0f;
  in.q_ref = 0.010f;
  in.tskin = 300.0f;
  in.q_surface = 0.012f;
  in.psfc = 95000.0f;

  const auto diag = surface::evaluate_neutral_surface_layer(in);
  TEST_CHECK(diag.wind_speed_ref > 0.0f);
  TEST_CHECK(diag.wind_speed_target > 0.0f);
  TEST_CHECK(diag.wind_speed_target < diag.wind_speed_ref);
  TEST_CHECK(diag.cm > 0.0f);
  TEST_CHECK(diag.ch > 0.0f);
  TEST_CHECK(diag.cq > 0.0f);
  TEST_CHECK(diag.ustar > 0.0f);
  TEST_CHECK(diag.t2 >= std::min(in.tskin, in.theta_ref));
  TEST_CHECK(diag.t2 <= std::max(in.tskin, in.theta_ref));
  TEST_CHECK(diag.q2 >= std::min(in.q_surface, in.q_ref));
  TEST_CHECK(diag.q2 <= std::max(in.q_surface, in.q_ref));
  TEST_CHECK(diag.rh2 >= 0.0f);
  TEST_CHECK(diag.rh2 <= 100.0f);

  auto rougher = in;
  rougher.z0m = 0.5f;
  const auto rougher_diag = surface::evaluate_neutral_surface_layer(rougher);
  TEST_CHECK(rougher_diag.wind_speed_target < diag.wind_speed_target);

  auto lower_target = in;
  lower_target.z_target_wind = 5.0f;
  const auto lower_target_diag =
      surface::evaluate_neutral_surface_layer(lower_target);
  TEST_CHECK(lower_target_diag.wind_speed_target < diag.wind_speed_target);

  return 0;
}
