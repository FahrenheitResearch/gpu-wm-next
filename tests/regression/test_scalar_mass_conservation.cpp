#include <algorithm>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/passive_tracer.hpp"
#include "gwm/dycore/dry_transport.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;
  auto layout = comm::VirtualRankLayout::build(32, 20, 1, 1, 2, 2, true, true);
  auto fields = comm::VirtualRankLayout::scatter_scalar(
      layout,
      [](int i, int j, int) {
        return static_cast<real>(1.0f + 0.1f * ((i + 2 * j) % 7));
      },
      "mass");

  double mass_before = 0.0;
  for (const auto& field : fields) {
    mass_before += field.owned_sum();
  }

  dycore::DryTransportConfig cfg{};
  cfg.u_adv = 12.0f;
  cfg.v_adv = 7.0f;
  cfg.dt = 2.0f;
  cfg.dx = 1000.0f;
  cfg.dy = 1000.0f;

  dycore::advance_scalar_ssprk3(fields, layout, cfg);

  double mass_after = 0.0;
  for (const auto& field : fields) {
    mass_after += field.owned_sum();
  }

  TEST_NEAR(mass_before, mass_after, 1.0e-3);

  domain::RectilinearDomainConfig domain_cfg{};
  domain_cfg.nx = 32;
  domain_cfg.ny = 20;
  domain_cfg.nz = 4;
  domain_cfg.halo = 1;
  domain_cfg.ranks_x = 2;
  domain_cfg.ranks_y = 2;
  domain_cfg.periodic_x = true;
  domain_cfg.periodic_y = true;
  const auto domain = domain::build_rectilinear_domain(domain_cfg);

  auto dry_states = dycore::make_constant_dry_state(domain.layout, 1.1f, 300.0f,
                                                    18.0f, -6.0f, 0.0f,
                                                    "tracer_dry");
  std::vector<real> qv_global(static_cast<std::size_t>(domain_cfg.nx) *
                                  static_cast<std::size_t>(domain_cfg.ny) *
                                  static_cast<std::size_t>(domain_cfg.nz),
                              0.0f);
  for (int k = 0; k < domain_cfg.nz; ++k) {
    for (int j = 0; j < domain_cfg.ny; ++j) {
      for (int i = 0; i < domain_cfg.nx; ++i) {
        const auto idx =
            (static_cast<std::size_t>(k) * static_cast<std::size_t>(domain_cfg.ny) +
             static_cast<std::size_t>(j)) *
                static_cast<std::size_t>(domain_cfg.nx) +
            static_cast<std::size_t>(i);
        qv_global[idx] = 0.003f + 0.0002f * static_cast<real>((i + 2 * j + k) % 9);
      }
    }
  }

  auto tracers = dycore::make_specific_humidity_tracers_from_global_field(
      dry_states, domain.layout, domain_cfg.nx, domain_cfg.ny, qv_global,
      "qv");

  double tracer_mass_before = 0.0;
  for (const auto& tracer_state : tracers) {
    tracer_mass_before += tracer_state.total_mass(0);
  }

  dycore::DryStepperConfig step_cfg{};
  step_cfg.dt = 2.0f;
  dycore::NullTracerBoundaryUpdater tracer_boundary;
  dycore::advance_passive_tracers_ssprk3(tracers, dry_states, domain.layout,
                                         domain.metrics, step_cfg,
                                         tracer_boundary);

  double tracer_mass_after = 0.0;
  real qv_min = 1.0e30f;
  real qv_max = -1.0e30f;
  for (std::size_t n = 0; n < tracers.size(); ++n) {
    tracer_mass_after += tracers[n].total_mass(0);
    const auto& rho_q =
        tracers[n].mass(gwm::state::kSpecificHumidityTracerName);
    const auto& rho_d = dry_states[n].rho_d;
    for (int k = 0; k < rho_q.nz(); ++k) {
      for (int j = 0; j < rho_q.ny(); ++j) {
        for (int i = 0; i < rho_q.nx(); ++i) {
          const real qv =
              rho_q(i, j, k) / std::max(rho_d(i, j, k), 1.0e-6f);
          qv_min = std::min(qv_min, qv);
          qv_max = std::max(qv_max, qv);
        }
      }
    }
  }

  TEST_NEAR(tracer_mass_before, tracer_mass_after, 1.0e-3);
  TEST_CHECK(qv_min >= 0.0f);
  TEST_CHECK(qv_max <= 0.0047f);
  return 0;
}
