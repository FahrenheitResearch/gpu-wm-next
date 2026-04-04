#include <algorithm>
#include <cmath>

#include "gwm/obsops/screen_obsops.hpp"
#include "gwm/surface/surface_layer_closure.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;
  obsops::ScreenInputs in{};
  in.z_ref = 40.0f;
  in.z_target_temp = 2.0f;
  in.z_target_wind = 10.0f;
  in.z0m = 0.1f;
  in.z0h = 0.02f;
  in.u_ref = 20.0f;
  in.theta_ref = 304.0f;
  in.q_ref = 0.010f;
  in.tskin = 300.0f;
  in.q_surface = 0.012f;
  in.psfc = 95000.0f;

  const auto diag = obsops::diagnose_neutral(in);
  const auto closure_diag = surface::evaluate_neutral_surface_layer(in);

  TEST_NEAR(diag.t2, closure_diag.t2, 1.0e-6f);
  TEST_NEAR(diag.q2, closure_diag.q2, 1.0e-6f);
  TEST_NEAR(diag.rh2, closure_diag.rh2, 1.0e-6f);
  TEST_NEAR(diag.u10, closure_diag.u10, 1.0e-6f);
  TEST_NEAR(diag.v10, closure_diag.v10, 1.0e-6f);
  TEST_NEAR(diag.wind_speed_ref, closure_diag.wind_speed_ref, 1.0e-6f);
  TEST_NEAR(diag.wind_speed_target, closure_diag.wind_speed_target, 1.0e-6f);
  TEST_NEAR(diag.cm, closure_diag.cm, 1.0e-6f);
  TEST_NEAR(diag.ch, closure_diag.ch, 1.0e-6f);
  TEST_NEAR(diag.cq, closure_diag.cq, 1.0e-6f);
  TEST_NEAR(diag.ustar, closure_diag.ustar, 1.0e-6f);
  TEST_CHECK(diag.u10 > 0.0f);
  TEST_CHECK(diag.u10 < in.u_ref);
  TEST_CHECK(diag.t2 >= std::min(in.tskin, in.theta_ref));
  TEST_CHECK(diag.t2 <= std::max(in.tskin, in.theta_ref));
  TEST_CHECK(diag.q2 >= std::min(in.q_surface, in.q_ref));
  TEST_CHECK(diag.q2 <= std::max(in.q_surface, in.q_ref));

  TEST_NEAR(obsops::saturation_specific_humidity(300.0f, in.psfc),
            surface::saturation_specific_humidity(300.0f, in.psfc), 1.0e-6f);

  auto moister = in;
  moister.q_ref = 0.014f;
  moister.q_surface = 0.016f;
  const auto moister_diag = obsops::diagnose_neutral(moister);
  TEST_CHECK(std::isfinite(moister_diag.q2));
  TEST_CHECK(std::isfinite(moister_diag.rh2));
  TEST_CHECK(moister_diag.q2 > diag.q2);
  TEST_CHECK(moister_diag.rh2 >= diag.rh2);
  TEST_CHECK(moister_diag.rh2 <= 100.0f);

  auto supersaturated = in;
  const auto qsat_surface =
      obsops::saturation_specific_humidity(in.tskin, in.psfc);
  supersaturated.q_ref = qsat_surface * 1.25f;
  supersaturated.q_surface = qsat_surface * 1.25f;
  const auto supersaturated_diag = obsops::diagnose_neutral(supersaturated);
  TEST_CHECK(supersaturated_diag.rh2 >= 0.0f);
  TEST_CHECK(supersaturated_diag.rh2 <= 100.0f);

  auto rougher = in;
  rougher.z0m = 0.5f;
  const auto rougher_diag = obsops::diagnose_neutral(rougher);
  TEST_CHECK(rougher_diag.u10 < diag.u10);
  return 0;
}
