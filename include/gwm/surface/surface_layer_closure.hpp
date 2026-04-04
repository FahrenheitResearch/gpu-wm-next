#pragma once

#include "gwm/core/types.hpp"

namespace gwm::surface {

struct SurfaceLayerInputs {
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

struct SurfaceLayerDiagnostics {
  real t2 = 0.0f;
  real q2 = 0.0f;
  real rh2 = 0.0f;
  real u10 = 0.0f;
  real v10 = 0.0f;
  real wind_speed_ref = 0.0f;
  real wind_speed_target = 0.0f;
  real cm = 0.0f;
  real ch = 0.0f;
  real cq = 0.0f;
  real ustar = 0.0f;
};

using SurfaceLayerOutputs = SurfaceLayerDiagnostics;

real saturation_specific_humidity(real temperature_k, real pressure_pa);
real relative_humidity_from_specific_humidity(real specific_humidity,
                                              real temperature_k,
                                              real pressure_pa);
real dewpoint_from_specific_humidity(real specific_humidity,
                                     real pressure_pa);

SurfaceLayerDiagnostics evaluate_neutral_surface_layer(
    const SurfaceLayerInputs& in);

}  // namespace gwm::surface
