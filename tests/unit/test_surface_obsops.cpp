#include "gwm/obsops/screen_obsops.hpp"

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
  TEST_CHECK(diag.u10 > 0.0f);
  TEST_CHECK(diag.u10 < in.u_ref);
  TEST_CHECK(diag.t2 >= std::min(in.tskin, in.theta_ref));
  TEST_CHECK(diag.t2 <= std::max(in.tskin, in.theta_ref));
  TEST_CHECK(diag.q2 >= std::min(in.q_surface, in.q_ref));
  TEST_CHECK(diag.q2 <= std::max(in.q_surface, in.q_ref));

  auto rougher = in;
  rougher.z0m = 0.5f;
  const auto rougher_diag = obsops::diagnose_neutral(rougher);
  TEST_CHECK(rougher_diag.u10 < diag.u10);
  return 0;
}
