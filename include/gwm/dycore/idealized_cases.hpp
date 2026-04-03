#pragma once

#include <string>
#include <vector>

#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"

namespace gwm::dycore {

struct ThermoBubbleConfig {
  real rho_background = 1.0f;
  real theta_background = 300.0f;
  real u_background = 0.0f;
  real v_background = 0.0f;
  real w_background = 0.0f;
  real center_x_fraction = 0.5f;
  real center_y_fraction = 0.5f;
  real center_z_fraction = 0.5f;
  real radius_x_fraction = 0.1f;
  real radius_y_fraction = 0.1f;
  real radius_z_fraction = 0.1f;
  real theta_perturbation = 2.0f;
};

struct MountainWaveBackgroundConfig {
  real rho_surface = 1.18f;
  real theta_ref = 300.0f;
  // Retained for compatibility with older configs; ignored by the current
  // constant-theta hydrostatic mountain-wave initializer.
  real density_scale_height = 8200.0f;
  real u_background = 15.0f;
  real v_background = 0.0f;
  real w_background = 0.0f;
};

std::vector<DryState> make_warm_bubble_state(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const ThermoBubbleConfig& config,
    const std::string& label_prefix);

std::vector<DryState> make_density_current_state(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const ThermoBubbleConfig& config,
    const std::string& label_prefix);

std::vector<DryState> make_mountain_wave_background_state(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics,
    const domain::IdealizedDomain& domain_topography,
    const MountainWaveBackgroundConfig& config,
    const std::string& label_prefix);

}  // namespace gwm::dycore
