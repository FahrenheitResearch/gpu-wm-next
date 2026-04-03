#include "gwm/dycore/idealized_cases.hpp"

#include <cmath>

#include "gwm/core/dry_thermo.hpp"

namespace gwm::dycore {

namespace {

constexpr real kPi = 3.14159265358979323846f;
constexpr real kGravity = 9.81f;

real bubble_weight(real dx, real dy, real dz) {
  const real r2 = dx * dx + dy * dy + dz * dz;
  if (r2 >= 1.0f) {
    return 0.0f;
  }
  return 0.5f * (1.0f + std::cos(kPi * std::sqrt(r2)));
}

real hydrostatic_rho_at_height(real rho_surface, real theta_ref, real z) {
  const real p_surface = gwm::core::dry_pressure_from_rho_theta_m(
      rho_surface, rho_surface * theta_ref);
  const real exner_surface =
      powf(fmaxf(p_surface / gwm::core::kReferencePressure, 1.0e-6f),
           gwm::core::kKappa);
  const real exner =
      fmaxf(1.0e-4f,
            exner_surface - kGravity * z / (gwm::core::dry_cp() * theta_ref));
  const real pressure =
      gwm::core::kReferencePressure * powf(exner, 1.0f / gwm::core::kKappa);
  return gwm::core::dry_rho_from_pressure_theta(pressure, theta_ref);
}

void set_background_face_momenta(DryState& state, real u_background,
                                 real v_background, real w_background) {
  auto face_density_x = [&](int i_face, int j, int k) {
    if (state.rho_d.nx() == 1) {
      return state.rho_d(0, j, k);
    }
    if (i_face <= 0) {
      return state.rho_d(0, j, k);
    }
    if (i_face >= state.rho_d.nx()) {
      return state.rho_d(state.rho_d.nx() - 1, j, k);
    }
    return 0.5f * (state.rho_d(i_face - 1, j, k) + state.rho_d(i_face, j, k));
  };

  auto face_density_y = [&](int i, int j_face, int k) {
    if (state.rho_d.ny() == 1) {
      return state.rho_d(i, 0, k);
    }
    if (j_face <= 0) {
      return state.rho_d(i, 0, k);
    }
    if (j_face >= state.rho_d.ny()) {
      return state.rho_d(i, state.rho_d.ny() - 1, k);
    }
    return 0.5f * (state.rho_d(i, j_face - 1, k) + state.rho_d(i, j_face, k));
  };

  auto face_density_z = [&](int i, int j, int k_face) {
    if (state.rho_d.nz() == 1) {
      return state.rho_d(i, j, 0);
    }
    if (k_face <= 0) {
      return state.rho_d(i, j, 0);
    }
    if (k_face >= state.rho_d.nz()) {
      return state.rho_d(i, j, state.rho_d.nz() - 1);
    }
    return 0.5f * (state.rho_d(i, j, k_face - 1) + state.rho_d(i, j, k_face));
  };

  auto& u_storage = state.mom_u.storage();
  for (int k = 0; k < u_storage.nz(); ++k) {
    for (int j = 0; j < u_storage.ny(); ++j) {
      for (int i = 0; i < u_storage.nx(); ++i) {
        u_storage(i, j, k) = face_density_x(i, j, k) * u_background;
      }
    }
  }

  auto& v_storage = state.mom_v.storage();
  for (int k = 0; k < v_storage.nz(); ++k) {
    for (int j = 0; j < v_storage.ny(); ++j) {
      for (int i = 0; i < v_storage.nx(); ++i) {
        v_storage(i, j, k) = face_density_y(i, j, k) * v_background;
      }
    }
  }

  auto& w_storage = state.mom_w.storage();
  for (int k = 0; k < w_storage.nz(); ++k) {
    for (int j = 0; j < w_storage.ny(); ++j) {
      for (int i = 0; i < w_storage.nx(); ++i) {
        w_storage(i, j, k) = face_density_z(i, j, k) * w_background;
      }
    }
  }
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
      for (int j = 0; j < desc.ny_local(); ++j) {
        const int j_global = desc.j_begin + j;
        const real y =
            static_cast<real>(j_global) * static_cast<real>(metrics.dy);
        for (int i = 0; i < desc.nx_local(); ++i) {
          const int i_global = desc.i_begin + i;
          const real x =
              static_cast<real>(i_global) * static_cast<real>(metrics.dx);
          const real z = metrics.z_center(i_global, j_global, k);
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
    const domain::IdealizedDomain& domain_topography,
    const MountainWaveBackgroundConfig& config,
    const std::string& label_prefix) {
  std::vector<DryState> states;
  states.reserve(layout.size());

  for (const auto& desc : layout) {
    DryState state(desc.nx_local(), desc.ny_local(), desc.nz, desc.halo,
                   label_prefix + "_rank_" + std::to_string(desc.rank));
    state.fill_constant(0.0f, config.theta_ref, config.u_background,
                        config.v_background, config.w_background);

    for (int k = 0; k < desc.nz; ++k) {
      for (int j = 0; j < desc.ny_local(); ++j) {
        const int j_global = desc.j_begin + j;
        for (int i = 0; i < desc.nx_local(); ++i) {
          const int i_global = desc.i_begin + i;
          const real terrain = domain_topography.terrain(i_global, j_global);
          const double z = metrics.flat_z_center(k) +
                           metrics.terrain_weight_center(k) *
                               static_cast<double>(terrain);
          (void)config.density_scale_height;
          const real rho_k =
              hydrostatic_rho_at_height(config.rho_surface, config.theta_ref,
                                        static_cast<real>(z));
          state.rho_d(i, j, k) = rho_k;
          state.rho_theta_m(i, j, k) = rho_k * config.theta_ref;
        }
      }
    }

    set_background_face_momenta(state, config.u_background, config.v_background,
                                config.w_background);

    states.push_back(std::move(state));
  }

  return states;
}

}  // namespace gwm::dycore
