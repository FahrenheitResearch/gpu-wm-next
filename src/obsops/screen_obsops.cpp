#include "gwm/obsops/screen_obsops.hpp"

namespace gwm::obsops {

ScreenDiagnostics diagnose_neutral(const ScreenInputs& in) {
  // Compatibility wrapper only; the shared closure owns the math.
  return surface::evaluate_neutral_surface_layer(in);
}

}  // namespace gwm::obsops
