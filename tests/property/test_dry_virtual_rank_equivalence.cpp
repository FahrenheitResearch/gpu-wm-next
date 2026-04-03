#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"

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

    double serial_mass = 0.0;
    double split_mass = 0.0;
    double serial_theta = 0.0;
    double split_theta = 0.0;
    for (const auto& s : serial_state) {
      serial_mass += s.total_dry_mass();
      serial_theta += s.total_rho_theta_m();
    }
    for (const auto& s : split_state) {
      split_mass += s.total_dry_mass();
      split_theta += s.total_rho_theta_m();
    }

    TEST_NEAR(serial_mass, split_mass, 1.0e-4);
    TEST_NEAR(serial_theta, split_theta, 1.0e-3);
    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
