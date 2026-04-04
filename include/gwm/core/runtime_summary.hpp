#pragma once

#include <string>
#include <vector>

#include "gwm/dycore/dry_diagnostics.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/physics/warm_rain.hpp"
#include "gwm/state/tracer_state.hpp"

namespace gwm::core {

struct TracerSummary {
  std::string name;
  std::string units;
  bool positive = true;
  double total_mass = 0.0;
  dycore::ScalarRange mixing_ratio;
};

struct MoistureBudgetSummary {
  double vapor_water_mass = 0.0;
  double cloud_water_mass = 0.0;
  double rain_water_mass = 0.0;
  double condensed_water_mass = 0.0;
  double total_water_mass = 0.0;
  double accumulated_surface_precipitation_sum_mm = 0.0;
  double mean_surface_precipitation_mm = 0.0;
  double max_surface_precipitation_mm = 0.0;
};

struct RuntimeStateSummary {
  dycore::DryStateSummary dry;
  std::vector<TracerSummary> tracers;
  MoistureBudgetSummary moisture;
};

[[nodiscard]] RuntimeStateSummary summarize_runtime_state(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::vector<physics::WarmRainSurfaceAccumulation>* accumulations =
        nullptr);

[[nodiscard]] std::string runtime_state_summary_to_json(
    const RuntimeStateSummary& summary, const std::string& indent = "  ");

}  // namespace gwm::core
