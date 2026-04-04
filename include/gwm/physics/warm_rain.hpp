#pragma once

#include <string>
#include <vector>

#include "gwm/domain/grid_metrics.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/state/field3d.hpp"
#include "gwm/state/tracer_state.hpp"

namespace gwm::physics {

struct WarmRainConfig {
  real dt = 2.0f;
  real condensation_relaxation = 1.0f;
  real cloud_autoconversion_threshold = 1.0e-4f;
  real cloud_autoconversion_rate = 1.0e-3f;
  real rain_evaporation_rate = 5.0e-3f;
  real rain_terminal_velocity = 0.0f;
  real sedimentation_layer_depth_m = 0.0f;
  bool enable_latent_heating = true;
};

struct WarmRainSurfaceAccumulation {
  int nx = 0;
  int ny = 0;
  double cell_area_m2 = 1.0;
  state::Field3D<real> liquid_precipitation_kg_m2;

  void reset(int nx_in, int ny_in, double cell_area_in = 1.0) {
    nx = nx_in;
    ny = ny_in;
    cell_area_m2 = cell_area_in;
    liquid_precipitation_kg_m2 =
        state::Field3D<real>(nx, ny, 1, 0, "warm_rain_surface_precip");
    liquid_precipitation_kg_m2.fill(0.0f);
  }

  [[nodiscard]] bool matches(int nx_in, int ny_in) const {
    return nx == nx_in && ny == ny_in && liquid_precipitation_kg_m2.nx() == nx &&
           liquid_precipitation_kg_m2.ny() == ny &&
           liquid_precipitation_kg_m2.nz() == 1 &&
           liquid_precipitation_kg_m2.halo() == 0;
  }
};

struct WarmRainSummary {
  double total_water_vapor_mass = 0.0;
  double total_cloud_water_mass = 0.0;
  double total_rain_water_mass = 0.0;
  double accumulated_surface_precipitation_sum_mm = 0.0;
  double mean_surface_precipitation_mm = 0.0;
  double max_surface_precipitation_mm = 0.0;
};

void apply_warm_rain_microphysics(std::vector<dycore::DryState>& states,
                                  std::vector<state::TracerState>& tracers,
                                  const WarmRainConfig& config,
                                  std::vector<WarmRainSurfaceAccumulation>*
                                      accumulations = nullptr);

void apply_warm_rain_microphysics(
    std::vector<dycore::DryState>& states,
    std::vector<state::TracerState>& tracers,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const WarmRainConfig& config,
    std::vector<WarmRainSurfaceAccumulation>* accumulations);

[[nodiscard]] WarmRainSummary summarize_warm_rain(
    const std::vector<state::TracerState>& tracers,
    const std::vector<WarmRainSurfaceAccumulation>* accumulations = nullptr);

[[nodiscard]] std::string warm_rain_summary_to_json(
    const WarmRainSummary& summary, const std::string& indent = "");

}  // namespace gwm::physics
