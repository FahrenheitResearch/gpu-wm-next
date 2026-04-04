#pragma once

#include <string>
#include <vector>

#include "gwm/domain/grid_metrics.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/dycore/dry_diagnostics.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/physics/warm_rain.hpp"
#include "gwm/state/tracer_state.hpp"

namespace gwm::core {

struct TracerSummary {
  std::string name;
  std::string units;
  bool positive = true;
  double column_integrated_sum_kg_m2 = 0.0;
  dycore::ScalarRange mixing_ratio;
};

struct MoistureBudgetSummary {
  double vapor_water_path_sum_kg_m2 = 0.0;
  double cloud_water_path_sum_kg_m2 = 0.0;
  double rain_water_path_sum_kg_m2 = 0.0;
  double condensed_water_path_sum_kg_m2 = 0.0;
  double total_water_path_sum_kg_m2 = 0.0;
  double accumulated_surface_precipitation_sum_mm = 0.0;
  double mean_surface_precipitation_mm = 0.0;
  double max_surface_precipitation_mm = 0.0;
  int total_surface_cell_count = 0;
  int precipitating_surface_cell_count = 0;
  double precipitating_surface_fraction = 0.0;
  double mean_precipitating_surface_precipitation_mm = 0.0;
};

struct RuntimeStateSummary {
  dycore::DryStateSummary dry;
  std::vector<TracerSummary> tracers;
  MoistureBudgetSummary moisture;
};

[[nodiscard]] RuntimeStateSummary summarize_runtime_state(
    const std::vector<dycore::DryState>& states,
    const std::vector<state::TracerState>& tracers,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics,
    const std::vector<physics::WarmRainSurfaceAccumulation>* accumulations =
        nullptr);

[[nodiscard]] std::string runtime_state_summary_to_json(
    const RuntimeStateSummary& summary, const std::string& indent = "  ");

}  // namespace gwm::core
