#pragma once

#include "gwm/surface/surface_layer_closure.hpp"

namespace gwm::obsops {

using ScreenInputs = surface::SurfaceLayerInputs;
using ScreenDiagnostics = surface::SurfaceLayerOutputs;

inline real saturation_specific_humidity(real temperature_k, real pressure_pa) {
  return surface::saturation_specific_humidity(temperature_k, pressure_pa);
}

// Screen obsops stay as a thin compatibility wrapper over the shared
// surface-layer closure so diagnostics and future flux work use one path.
ScreenDiagnostics diagnose_neutral(const ScreenInputs& in);

}  // namespace gwm::obsops
