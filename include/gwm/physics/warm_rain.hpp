#pragma once

#include <string>
#include <vector>

#include "gwm/dycore/dry_core.hpp"
#include "gwm/state/tracer_state.hpp"

namespace gwm::physics {

struct WarmRainConfig {
  real dt = 2.0f;
  real condensation_relaxation = 1.0f;
  real cloud_autoconversion_threshold = 1.0e-4f;
  real cloud_autoconversion_rate = 1.0e-3f;
  real rain_evaporation_rate = 5.0e-3f;
  bool enable_latent_heating = true;
};

struct WarmRainSummary {
  double total_water_vapor_mass = 0.0;
  double total_cloud_water_mass = 0.0;
  double total_rain_water_mass = 0.0;
};

void apply_warm_rain_microphysics(std::vector<dycore::DryState>& states,
                                  std::vector<state::TracerState>& tracers,
                                  const WarmRainConfig& config);

[[nodiscard]] WarmRainSummary summarize_warm_rain(
    const std::vector<state::TracerState>& tracers);

[[nodiscard]] std::string warm_rain_summary_to_json(
    const WarmRainSummary& summary, const std::string& indent = "");

}  // namespace gwm::physics
