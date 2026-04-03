#include <cmath>

#include "gwm/comm/halo_exchange.hpp"
#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/dry_momentum_flux.hpp"
#include "gwm/dycore/dry_pressure_gradient.hpp"
#include "gwm/dycore/idealized_cases.hpp"

#include "test_assert.hpp"

namespace {

using gwm::real;

std::vector<gwm::state::Field3D<real>> clone_rho_fields(
    const std::vector<gwm::dycore::DryState>& states) {
  std::vector<gwm::state::Field3D<real>> fields;
  fields.reserve(states.size());
  for (const auto& state : states) {
    fields.push_back(state.rho_d.clone_empty_like("_rho"));
    fields.back().copy_all_from(state.rho_d);
  }
  return fields;
}

std::vector<gwm::state::Field3D<real>> clone_theta_fields(
    const std::vector<gwm::dycore::DryState>& states) {
  std::vector<gwm::state::Field3D<real>> fields;
  fields.reserve(states.size());
  for (const auto& state : states) {
    fields.push_back(state.rho_theta_m.clone_empty_like("_theta"));
    fields.back().copy_all_from(state.rho_theta_m);
  }
  return fields;
}

std::vector<gwm::state::FaceField<real>> clone_face_fields(
    const std::vector<gwm::dycore::DryState>& states,
    gwm::state::FaceOrientation orientation) {
  std::vector<gwm::state::FaceField<real>> fields;
  fields.reserve(states.size());
  for (const auto& state : states) {
    const auto& src = orientation == gwm::state::FaceOrientation::X
                          ? state.mom_u
                          : orientation == gwm::state::FaceOrientation::Y
                                ? state.mom_v
                                : state.mom_w;
    fields.emplace_back(state.rho_d.nx(), state.rho_d.ny(), state.rho_d.nz(),
                        state.rho_d.halo(), orientation, "face_clone");
    fields.back().storage().copy_all_from(src.storage());
  }
  return fields;
}

std::vector<gwm::dycore::DryState> clone_states_with_zero_momentum(
    const std::vector<gwm::dycore::DryState>& states) {
  std::vector<gwm::dycore::DryState> out;
  out.reserve(states.size());
  for (const auto& state : states) {
    out.push_back(state.clone_empty_like("zero_mom"));
    out.back().copy_all_from(state);
    out.back().mom_u.storage().fill(0.0f);
    out.back().mom_v.storage().fill(0.0f);
    out.back().mom_w.storage().fill(0.0f);
  }
  return out;
}

std::vector<gwm::dycore::DrySlowTendencies> make_tendency_bundle(
    const std::vector<gwm::dycore::DryState>& states, real u_fill = 0.0f,
    real v_fill = 0.0f, real w_fill = 0.0f) {
  std::vector<gwm::dycore::DrySlowTendencies> out;
  out.reserve(states.size());
  for (const auto& state : states) {
    out.emplace_back(state.rho_d.nx(), state.rho_d.ny(), state.rho_d.nz(),
                     state.rho_d.halo(), "tendency");
    out.back().fill_zero();
    out.back().mom_u.storage().fill(u_fill);
    out.back().mom_v.storage().fill(v_fill);
    out.back().mom_w.storage().fill(w_fill);
  }
  return out;
}

void check_face_offset_preserved(const gwm::state::FaceField<real>& shifted,
                                 const gwm::state::FaceField<real>& base,
                                 real expected, real tol) {
  const auto& a = shifted.storage();
  const auto& b = base.storage();
  for (int k = 0; k < a.nz(); ++k) {
    for (int j = 0; j < a.ny(); ++j) {
      for (int i = 0; i < a.nx(); ++i) {
        TEST_NEAR(a(i, j, k) - b(i, j, k), expected, tol);
      }
    }
  }
}

double max_face_abs_diff(const gwm::state::FaceField<real>& a,
                         const gwm::state::FaceField<real>& b) {
  double max_diff = 0.0;
  const auto& sa = a.storage();
  const auto& sb = b.storage();
  for (int k = 0; k < sa.nz(); ++k) {
    for (int j = 0; j < sa.ny(); ++j) {
      for (int i = 0; i < sa.nx(); ++i) {
        max_diff = std::max(max_diff,
                            std::fabs(static_cast<double>(sa(i, j, k) - sb(i, j, k))));
      }
    }
  }
  return max_diff;
}

}  // namespace

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig cfg{};
    cfg.nx = 24;
    cfg.ny = 20;
    cfg.nz = 10;
    cfg.halo = 1;
    cfg.ranks_x = 2;
    cfg.ranks_y = 2;
    cfg.terrain_kind = domain::TerrainProfileKind::CosineMountain;
    cfg.mountain_height = 900.0;
    cfg.mountain_half_width_x = 0.16;
    cfg.mountain_half_width_y = 0.28;

    const auto domain = domain::build_rectilinear_domain(cfg);

    dycore::ThermoBubbleConfig bubble_cfg{};
    bubble_cfg.rho_background = 1.0f;
    bubble_cfg.theta_background = 300.0f;
    bubble_cfg.u_background = 8.0f;
    bubble_cfg.v_background = -3.0f;
    bubble_cfg.w_background = 0.35f;
    bubble_cfg.theta_perturbation = 3.0f;
    bubble_cfg.center_x_fraction = 0.5f;
    bubble_cfg.center_y_fraction = 0.5f;
    bubble_cfg.center_z_fraction = 0.35f;
    bubble_cfg.radius_x_fraction = 0.15f;
    bubble_cfg.radius_y_fraction = 0.15f;
    bubble_cfg.radius_z_fraction = 0.18f;

    auto composed_states = dycore::make_warm_bubble_state(
        domain.layout, domain.metrics, bubble_cfg, "composition");
    for (auto& state : composed_states) {
      auto& u = state.mom_u.storage();
      for (int k = 0; k < u.nz(); ++k) {
        for (int j = 0; j < u.ny(); ++j) {
          for (int i = 0; i < u.nx(); ++i) {
            const real phase = static_cast<real>(
                0.05 * std::sin(0.35 * static_cast<double>(i)) +
                0.03 * std::cos(0.40 * static_cast<double>(j)));
            u(i, j, k) *= (1.0f + phase);
          }
        }
      }

      auto& v = state.mom_v.storage();
      for (int k = 0; k < v.nz(); ++k) {
        for (int j = 0; j < v.ny(); ++j) {
          for (int i = 0; i < v.nx(); ++i) {
            const real phase = static_cast<real>(
                0.04 * std::cos(0.25 * static_cast<double>(i)) -
                0.02 * std::sin(0.30 * static_cast<double>(j)));
            v(i, j, k) *= (1.0f + phase);
          }
        }
      }

      auto& w = state.mom_w.storage();
      for (int k = 0; k < w.nz(); ++k) {
        for (int j = 0; j < w.ny(); ++j) {
          for (int i = 0; i < w.nx(); ++i) {
            const real phase = static_cast<real>(
                0.03 * std::sin(0.20 * static_cast<double>(i)) +
                0.02 * std::cos(0.27 * static_cast<double>(j)) +
                0.01 * std::sin(0.31 * static_cast<double>(k)));
            w(i, j, k) *= (1.0f + phase);
          }
        }
      }
    }

    auto rho_fields = clone_rho_fields(composed_states);
    auto theta_fields = clone_theta_fields(composed_states);
    auto mom_u_fields =
        clone_face_fields(composed_states, state::FaceOrientation::X);
    auto mom_v_fields =
        clone_face_fields(composed_states, state::FaceOrientation::Y);
    auto mom_w_fields =
        clone_face_fields(composed_states, state::FaceOrientation::Z);

    comm::HaloExchange::exchange_scalar(rho_fields, domain.layout);
    comm::HaloExchange::exchange_scalar(theta_fields, domain.layout);
    comm::HaloExchange::exchange_face(mom_u_fields, domain.layout);
    comm::HaloExchange::exchange_face(mom_v_fields, domain.layout);
    comm::HaloExchange::exchange_face(mom_w_fields, domain.layout);

    std::vector<state::Field3D<real>> pressure_fields;
    dycore::compute_dry_pressure_fields(rho_fields, theta_fields, pressure_fields);

    {
      auto base = make_tendency_bundle(composed_states);
      auto shifted = make_tendency_bundle(composed_states, 1.25f, -0.75f, 0.5f);

      dycore::add_dry_momentum_flux_tendencies(
          rho_fields, mom_u_fields, mom_v_fields, mom_w_fields, domain.layout,
          domain.metrics, base);
      dycore::add_dry_momentum_flux_tendencies(
          rho_fields, mom_u_fields, mom_v_fields, mom_w_fields, domain.layout,
          domain.metrics, shifted);
      dycore::add_horizontal_pressure_gradient_tendencies(
          pressure_fields, domain.layout, domain.metrics, base);
      dycore::add_horizontal_pressure_gradient_tendencies(
          pressure_fields, domain.layout, domain.metrics, shifted);

      for (std::size_t n = 0; n < base.size(); ++n) {
        check_face_offset_preserved(shifted[n].mom_u, base[n].mom_u, 1.25f,
                                    1.0e-5f);
        check_face_offset_preserved(shifted[n].mom_v, base[n].mom_v, -0.75f,
                                    1.0e-5f);
        check_face_offset_preserved(shifted[n].mom_w, base[n].mom_w, 0.5f,
                                    1.0e-5f);
      }
    }

    {
      auto pg_only = make_tendency_bundle(composed_states);
      auto full = make_tendency_bundle(composed_states);

      dycore::add_horizontal_pressure_gradient_tendencies(
          pressure_fields, domain.layout, domain.metrics, pg_only);
      dycore::add_dry_momentum_flux_tendencies(
          rho_fields, mom_u_fields, mom_v_fields, mom_w_fields, domain.layout,
          domain.metrics, full);
      dycore::add_horizontal_pressure_gradient_tendencies(
          pressure_fields, domain.layout, domain.metrics, full);

      double max_u_diff = 0.0;
      double max_v_diff = 0.0;
      for (std::size_t n = 0; n < full.size(); ++n) {
        max_u_diff = std::max(max_u_diff,
                              max_face_abs_diff(full[n].mom_u, pg_only[n].mom_u));
        max_v_diff = std::max(max_v_diff,
                              max_face_abs_diff(full[n].mom_v, pg_only[n].mom_v));
      }

      TEST_CHECK(max_u_diff > 1.0e-6);
      TEST_CHECK(max_v_diff > 1.0e-6);
    }

    {
      auto zero_mom_states = clone_states_with_zero_momentum(composed_states);
      std::vector<dycore::DrySlowTendencies> full_slow;
      std::vector<dycore::DrySlowTendencies> zero_mom_slow;
      dycore::DryStepperConfig step_cfg{};

      dycore::compute_slow_tendencies(composed_states, domain.layout,
                                      domain.metrics, step_cfg, full_slow);
      dycore::compute_slow_tendencies(zero_mom_states, domain.layout,
                                      domain.metrics, step_cfg, zero_mom_slow);

      double max_u_diff = 0.0;
      double max_v_diff = 0.0;
      double max_w_diff = 0.0;
      for (std::size_t n = 0; n < full_slow.size(); ++n) {
        max_u_diff =
            std::max(max_u_diff,
                     max_face_abs_diff(full_slow[n].mom_u, zero_mom_slow[n].mom_u));
        max_v_diff =
            std::max(max_v_diff,
                     max_face_abs_diff(full_slow[n].mom_v, zero_mom_slow[n].mom_v));
        max_w_diff =
            std::max(max_w_diff,
                     max_face_abs_diff(full_slow[n].mom_w, zero_mom_slow[n].mom_w));
      }

      TEST_CHECK(max_u_diff > 1.0e-6);
      TEST_CHECK(max_v_diff > 1.0e-6);
      TEST_CHECK(max_w_diff > 1.0e-6);
    }

    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
