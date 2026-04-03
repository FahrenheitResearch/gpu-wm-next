#include <cmath>
#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/dry_diagnostics.hpp"
#include "gwm/dycore/idealized_cases.hpp"

#include "test_assert.hpp"

int main() {
  try {
    using namespace gwm;

    domain::RectilinearDomainConfig flat{};
    flat.nx = 96;
    flat.ny = 24;
    flat.nz = 24;
    flat.halo = 1;
    flat.ranks_x = 2;
    flat.ranks_y = 2;
    flat.dx = 1000.0;
    flat.dy = 1000.0;
    flat.z_top = 18000.0;

    domain::RectilinearDomainConfig mountain = flat;
    mountain.terrain_kind = domain::TerrainProfileKind::CosineMountain;
    mountain.mountain_height = 1200.0;
    mountain.mountain_half_width_x = 0.10;
    mountain.mountain_half_width_y = 0.40;

    const auto flat_domain = domain::build_rectilinear_domain(flat);
    const auto mountain_domain = domain::build_rectilinear_domain(mountain);

    dycore::MountainWaveBackgroundConfig bg{};
    bg.rho_surface = 1.18f;
    bg.theta_ref = 300.0f;
    bg.density_scale_height = 8200.0f;
    bg.u_background = 15.0f;

    auto flat_states = dycore::make_mountain_wave_background_state(
        flat_domain.layout, flat_domain.metrics, flat_domain, bg, "flat_wave");
    auto mountain_states = dycore::make_mountain_wave_background_state(
        mountain_domain.layout, mountain_domain.metrics, mountain_domain, bg,
        "terrain_wave");

    double flat_mass_before = 0.0;
    double terrain_mass_before = 0.0;
    for (const auto& state : flat_states) {
      flat_mass_before += state.total_dry_mass();
    }
    for (const auto& state : mountain_states) {
      terrain_mass_before += state.total_dry_mass();
    }

    dycore::DryStepperConfig step_cfg{};
    step_cfg.dt = 1.0f;
    step_cfg.fast_substeps = 4;
    dycore::NullBoundaryUpdater boundary;
    dycore::LocalSplitExplicitFastMode fast_modes;

    for (int step = 0; step < 6; ++step) {
      dycore::advance_dry_state_ssprk3(flat_states, flat_domain.layout,
                                       flat_domain.metrics, step_cfg, boundary,
                                       fast_modes);
      dycore::advance_dry_state_ssprk3(mountain_states, mountain_domain.layout,
                                       mountain_domain.metrics, step_cfg,
                                       boundary, fast_modes);
    }

    double flat_mass_after = 0.0;
    double terrain_mass_after = 0.0;
    for (const auto& state : flat_states) {
      flat_mass_after += state.total_dry_mass();
    }
    for (const auto& state : mountain_states) {
      terrain_mass_after += state.total_dry_mass();
    }

    const auto flat_summary = dycore::summarize_dry_states(flat_states);
    const auto terrain_summary = dycore::summarize_dry_states(mountain_states);
    TEST_NEAR(flat_mass_before, flat_mass_after, 5.0e-2);
    TEST_NEAR(terrain_mass_before, terrain_mass_after, 1.0);
    TEST_CHECK(flat_summary.w_face.max < 2.0e-2);
    TEST_CHECK(flat_summary.w_face.min > -2.0e-2);
    TEST_CHECK(terrain_summary.w_face.max > 1.0e-2 ||
               terrain_summary.w_face.min < -1.0e-2);
    TEST_CHECK(std::isfinite(terrain_summary.total_dry_mass));
    return 0;
  } catch (const std::exception& ex) {
    test_fail(ex.what());
  }
}
