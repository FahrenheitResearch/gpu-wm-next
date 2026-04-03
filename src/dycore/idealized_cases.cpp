#include "gwm/dycore/idealized_cases.hpp"

#include <cmath>

namespace gwm::dycore {

namespace {

constexpr real kPi = 3.14159265358979323846f;

real bubble_weight(real dx, real dy, real dz) {
  const real r2 = dx * dx + dy * dy + dz * dz;
  if (r2 >= 1.0f) {
    return 0.0f;
  }
  return 0.5f * (1.0f + std::cos(kPi * std::sqrt(r2)));
}

std::vector<DryState> make_bubble_like_state(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const ThermoBubbleConfig& config,
    const std::string& label_prefix, real signed_theta_perturbation) {
  std::vector<DryState> states;
  states.reserve(layout.size());

  const real x0 = config.center_x_fraction * static_cast<real>(metrics.nx - 1) *
                  static_cast<real>(metrics.dx);
  const real y0 = config.center_y_fraction * static_cast<real>(metrics.ny - 1) *
                  static_cast<real>(metrics.dy);
  const real z0 = config.center_z_fraction * static_cast<real>(metrics.z_top);

  const real rx = std::max(config.radius_x_fraction *
                               static_cast<real>(metrics.nx) *
                               static_cast<real>(metrics.dx),
                           1.0f);
  const real ry = std::max(config.radius_y_fraction *
                               static_cast<real>(metrics.ny) *
                               static_cast<real>(metrics.dy),
                           1.0f);
  const real rz = std::max(config.radius_z_fraction *
                               static_cast<real>(metrics.z_top),
                           1.0f);

  for (const auto& desc : layout) {
    DryState state(desc.nx_local(), desc.ny_local(), desc.nz, desc.halo,
                   label_prefix + "_rank_" + std::to_string(desc.rank));
    state.fill_constant(config.rho_background, config.theta_background,
                        config.u_background, config.v_background,
                        config.w_background);

    for (int k = 0; k < desc.nz; ++k) {
      const real z = static_cast<real>(metrics.z_centers[static_cast<std::size_t>(k)]);
      for (int j = 0; j < desc.ny_local(); ++j) {
        const real y =
            static_cast<real>(desc.j_begin + j) * static_cast<real>(metrics.dy);
        for (int i = 0; i < desc.nx_local(); ++i) {
          const real x = static_cast<real>(desc.i_begin + i) *
                         static_cast<real>(metrics.dx);
          const real wx = (x - x0) / rx;
          const real wy = (y - y0) / ry;
          const real wz = (z - z0) / rz;
          const real weight = bubble_weight(wx, wy, wz);
          const real theta =
              config.theta_background + signed_theta_perturbation * weight;
          state.rho_theta_m(i, j, k) = config.rho_background * theta;
        }
      }
    }

    states.push_back(std::move(state));
  }

  return states;
}

}  // namespace

std::vector<DryState> make_warm_bubble_state(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const ThermoBubbleConfig& config,
    const std::string& label_prefix) {
  return make_bubble_like_state(layout, metrics, config, label_prefix,
                                std::abs(config.theta_perturbation));
}

std::vector<DryState> make_density_current_state(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const ThermoBubbleConfig& config,
    const std::string& label_prefix) {
  return make_bubble_like_state(layout, metrics, config, label_prefix,
                                -std::abs(config.theta_perturbation));
}

std::vector<DryState> make_mountain_wave_background_state(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics,
    const domain::IdealizedDomain& /*domain_topography*/,
    const MountainWaveBackgroundConfig& config,
    const std::string& label_prefix) {
  auto states = make_hydrostatic_rest_state(layout, metrics, config.rho_surface,
                                            config.theta_ref,
                                            config.density_scale_height,
                                            label_prefix);

  for (auto& state : states) {
    state.mom_u.storage().fill(config.rho_surface * config.u_background);
    state.mom_v.storage().fill(config.rho_surface * config.v_background);
    state.mom_w.storage().fill(config.rho_surface * config.w_background);
  }

  return states;
}

}  // namespace gwm::dycore
