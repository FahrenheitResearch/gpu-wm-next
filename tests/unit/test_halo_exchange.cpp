#include "gwm/comm/halo_exchange.hpp"
#include "gwm/comm/virtual_rank_layout.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;
  auto layout = comm::VirtualRankLayout::build(8, 8, 1, 1, 2, 2, true, true);
  auto fields = comm::VirtualRankLayout::scatter_scalar(
      layout, [](int i, int j, int) { return static_cast<real>(100 * j + i); },
      "halo");

  comm::HaloExchange::exchange_scalar(fields, layout);

  const auto& rank00 = fields[0];
  TEST_NEAR(rank00(-1, 0, 0), 7.0f, 1.0e-6f);
  TEST_NEAR(rank00(4, 0, 0), 4.0f, 1.0e-6f);
  TEST_NEAR(rank00(0, -1, 0), 700.0f, 1.0e-6f);
  TEST_NEAR(rank00(0, 4, 0), 400.0f, 1.0e-6f);
  return 0;
}
