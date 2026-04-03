#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/dycore/dry_transport.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;
  auto layout = comm::VirtualRankLayout::build(24, 24, 3, 1, 2, 2, true, true);
  auto fields = comm::VirtualRankLayout::scatter_scalar(
      layout, [](int, int, int) { return 7.5f; }, "constant");

  dycore::DryTransportConfig cfg{};
  cfg.u_adv = 20.0f;
  cfg.v_adv = -10.0f;
  cfg.dt = 4.0f;
  cfg.dx = 1000.0f;
  cfg.dy = 1000.0f;

  dycore::advance_scalar_ssprk3(fields, layout, cfg);
  const auto gathered = comm::VirtualRankLayout::gather_scalar(fields, layout);
  for (int k = 0; k < gathered.nz(); ++k) {
    for (int j = 0; j < gathered.ny(); ++j) {
      for (int i = 0; i < gathered.nx(); ++i) {
        TEST_NEAR(gathered(i, j, k), 7.5f, 1.0e-5f);
      }
    }
  }
  return 0;
}
