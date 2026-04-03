#include "gwm/comm/virtual_rank_layout.hpp"
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
  return 0;
}
