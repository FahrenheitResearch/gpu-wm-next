#pragma once

#include <algorithm>
#include <cmath>

#include "gwm/core/types.hpp"

namespace gwm::obsops {

struct ScreenInputs {
  real z_ref = 40.0f;
  real z_target_temp = 2.0f;
  real z_target_wind = 10.0f;
  real z0m = 0.1f;
  real z0h = 0.02f;
  real u_ref = 0.0f;
  real v_ref = 0.0f;
  real theta_ref = 300.0f;
  real q_ref = 0.010f;
  real tskin = 299.0f;
  real q_surface = 0.011f;
  real psfc = 95000.0f;
};

struct ScreenDiagnostics {
  real t2 = 0.0f;
  real q2 = 0.0f;
  real rh2 = 0.0f;
  real u10 = 0.0f;
  real v10 = 0.0f;
};

real saturation_specific_humidity(real temperature_k, real pressure_pa);
ScreenDiagnostics diagnose_neutral(const ScreenInputs& in);

}  // namespace gwm::obsops
