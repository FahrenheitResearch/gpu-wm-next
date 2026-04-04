#pragma once

#include <vector>

#include "gwm/surface/surface_layer_closure.hpp"
#include "gwm/surface/surface_state.hpp"
#include "gwm/surface/surface_static_properties.hpp"

namespace gwm::surface {

struct SurfaceLayerExchangeForcing {
  real z_ref = 40.0f;
  real z_target_temp = 2.0f;
  real z_target_wind = 10.0f;
  real u_ref = 0.0f;
  real v_ref = 0.0f;
  real theta_ref = 300.0f;
  real q_ref = 0.010f;
  real psfc = 95000.0f;
  std::vector<real> tile_q_surface;
};

struct SurfaceLayerExchangeResult {
  SurfaceLayerDiagnostics cell_mean{};
  std::vector<SurfaceLayerDiagnostics> tile_outputs;
};

SurfaceLayerExchangeResult evaluate_neutral_surface_exchange(
    const SurfaceState& state, const SurfaceStaticProperties& properties, int i,
    int j, const SurfaceLayerExchangeForcing& forcing);

}  // namespace gwm::surface
