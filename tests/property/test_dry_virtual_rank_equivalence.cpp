#include <algorithm>
#include <cmath>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/passive_tracer.hpp"
#include "gwm/dycore/dry_transport.hpp"

#include "test_assert.hpp"

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig serial_cfg{};
    serial_cfg.nx = 24;
    serial_cfg.ny = 20;
    serial_cfg.nz = 8;
    serial_cfg.halo = 1;
    serial_cfg.ranks_x = 1;
    serial_cfg.ranks_y = 1;

    domain::RectilinearDomainConfig split_cfg = serial_cfg;
    split_cfg.ranks_x = 2;
    split_cfg.ranks_y = 2;

    const auto serial_domain = domain::build_rectilinear_domain(serial_cfg);
    const auto split_domain = domain::build_rectilinear_domain(split_cfg);

    auto serial_state = dycore::make_constant_dry_state(
        serial_domain.layout, 1.15f, 302.0f, 12.0f, -6.0f, 0.0f, "serial");
    auto split_state = dycore::make_constant_dry_state(
        split_domain.layout, 1.15f, 302.0f, 12.0f, -6.0f, 0.0f, "split");

    dycore::DryStepperConfig cfg{};
    cfg.dt = 2.0f;
    cfg.fast_substeps = 2;
    dycore::NullBoundaryUpdater boundary;
    dycore::LocalSplitExplicitFastMode fast_modes;

    dycore::advance_dry_state_ssprk3(serial_state, serial_domain.layout,
                                     serial_domain.metrics, cfg, boundary,
                                     fast_modes);
    dycore::advance_dry_state_ssprk3(split_state, split_domain.layout,
                                     split_domain.metrics, cfg, boundary,
                                     fast_modes);

    std::vector<real> qv_global(static_cast<std::size_t>(serial_cfg.nx) *
                                    static_cast<std::size_t>(serial_cfg.ny) *
                                    static_cast<std::size_t>(serial_cfg.nz),
                                0.0f);
    for (int k = 0; k < serial_cfg.nz; ++k) {
      for (int j = 0; j < serial_cfg.ny; ++j) {
        for (int i = 0; i < serial_cfg.nx; ++i) {
          const auto idx =
              (static_cast<std::size_t>(k) *
                   static_cast<std::size_t>(serial_cfg.ny) +
               static_cast<std::size_t>(j)) *
                  static_cast<std::size_t>(serial_cfg.nx) +
              static_cast<std::size_t>(i);
          qv_global[idx] =
              0.0025f + 0.0001f * static_cast<real>((i + j + 2 * k) % 11);
        }
      }
    }

    auto serial_tracers = dycore::make_specific_humidity_tracers_from_global_field(
        serial_state, serial_domain.layout, serial_cfg.nx, serial_cfg.ny,
        qv_global, "serial_qv");
    auto split_tracers = dycore::make_specific_humidity_tracers_from_global_field(
        split_state, split_domain.layout, split_cfg.nx, split_cfg.ny, qv_global,
        "split_qv");

    dycore::NullTracerBoundaryUpdater tracer_boundary;
    dycore::advance_passive_tracers_ssprk3(serial_tracers, serial_state,
                                           serial_domain.layout,
                                           serial_domain.metrics, cfg,
                                           tracer_boundary);
    dycore::advance_passive_tracers_ssprk3(split_tracers, split_state,
                                           split_domain.layout,
                                           split_domain.metrics, cfg,
                                           tracer_boundary);

    double serial_mass = 0.0;
    double split_mass = 0.0;
    double serial_theta = 0.0;
    double split_theta = 0.0;
    double serial_tracer_mass = 0.0;
    double split_tracer_mass = 0.0;
    for (const auto& s : serial_state) {
      serial_mass += s.total_dry_mass();
      serial_theta += s.total_rho_theta_m();
    }
    for (const auto& s : split_state) {
      split_mass += s.total_dry_mass();
      split_theta += s.total_rho_theta_m();
    }
    for (const auto& t : serial_tracers) {
      serial_tracer_mass += t.total_mass(0);
    }
    for (const auto& t : split_tracers) {
      split_tracer_mass += t.total_mass(0);
    }

    TEST_NEAR(serial_mass, split_mass, 1.0e-4);
    TEST_NEAR(serial_theta, split_theta, 1.0e-3);

    std::vector<state::Field3D<real>> split_qv_fields;
    split_qv_fields.reserve(split_tracers.size());
    for (const auto& tracer_state : split_tracers) {
      split_qv_fields.push_back(tracer_state.mass(0).clone_empty_like("_gather"));
      split_qv_fields.back().copy_all_from(tracer_state.mass(0));
    }
    const auto gathered_split_qv =
        comm::VirtualRankLayout::gather_scalar(split_qv_fields,
                                               split_domain.layout,
                                               "split_qv_global");
    const auto& serial_qv = serial_tracers.front().mass(0);
    real max_qv_diff = 0.0f;
    for (int k = 0; k < serial_qv.nz(); ++k) {
      for (int j = 0; j < serial_qv.ny(); ++j) {
        for (int i = 0; i < serial_qv.nx(); ++i) {
          max_qv_diff = std::max(
              max_qv_diff,
              std::fabs(serial_qv(i, j, k) - gathered_split_qv(i, j, k)));
        }
      }
    }

    TEST_NEAR(serial_tracer_mass, split_tracer_mass, 1.0e-3);
    TEST_CHECK(max_qv_diff < 1.0e-4f);
    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
