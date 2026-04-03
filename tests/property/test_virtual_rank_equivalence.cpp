#include <cmath>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/dycore/dry_transport.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;
  auto init = [](int i, int j, int) {
    return static_cast<real>(std::sin(0.2 * i) + std::cos(0.3 * j));
  };

  auto layout_serial = comm::VirtualRankLayout::build(32, 32, 2, 1, 1, 1, true, true);
  auto fields_serial =
      comm::VirtualRankLayout::scatter_scalar(layout_serial, init, "serial");

  auto layout_split = comm::VirtualRankLayout::build(32, 32, 2, 1, 2, 2, true, true);
  auto fields_split =
      comm::VirtualRankLayout::scatter_scalar(layout_split, init, "split");

  dycore::DryTransportConfig cfg{};
  cfg.u_adv = 15.0f;
  cfg.v_adv = -5.0f;
  cfg.dt = 5.0f;
  cfg.dx = 1000.0f;
  cfg.dy = 1000.0f;

  dycore::advance_scalar_ssprk3(fields_serial, layout_serial, cfg);
  dycore::advance_scalar_ssprk3(fields_split, layout_split, cfg);

  const auto serial = comm::VirtualRankLayout::gather_scalar(fields_serial, layout_serial);
  const auto split = comm::VirtualRankLayout::gather_scalar(fields_split, layout_split);

  double max_abs_err = 0.0;
  for (int k = 0; k < serial.nz(); ++k) {
    for (int j = 0; j < serial.ny(); ++j) {
      for (int i = 0; i < serial.nx(); ++i) {
        max_abs_err = std::max(
            max_abs_err,
            std::fabs(static_cast<double>(serial(i, j, k) - split(i, j, k))));
      }
    }
  }

  TEST_NEAR(max_abs_err, 0.0, 1.0e-5);
  return 0;
}
