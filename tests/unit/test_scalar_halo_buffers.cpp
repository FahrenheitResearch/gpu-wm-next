#include "gwm/comm/cartesian_topology.hpp"
#include "gwm/comm/halo_exchange.hpp"
#include "gwm/comm/virtual_rank_layout.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;

  auto layout = comm::VirtualRankLayout::build(8, 8, 2, 1, 2, 2, true, true);
  auto fields = comm::VirtualRankLayout::scatter_scalar(
      layout, [](int i, int j, int k) { return static_cast<real>(1000 * k + 100 * j + i); },
      "buffers");

  const auto plan = comm::make_scalar_halo_exchange_plan(layout, layout[0]);
  TEST_CHECK(plan.neighbors.west == 1);
  TEST_CHECK(plan.neighbors.east == 1);
  TEST_CHECK(plan.neighbors.south == 2);
  TEST_CHECK(plan.neighbors.north == 2);
  TEST_CHECK(plan.x_face_count == 8);
  TEST_CHECK(plan.y_face_count == 12);

  auto buffers = comm::HaloExchange::allocate_scalar_buffers(plan);
  comm::HaloExchange::pack_scalar_x_faces(fields[0], layout[0], buffers);
  TEST_NEAR(buffers.west_send[0], 0.0f, 1.0e-6f);
  TEST_NEAR(buffers.east_send[0], 3.0f, 1.0e-6f);

  comm::HaloExchange::exchange_scalar(fields, layout);
  TEST_NEAR(fields[0](-1, 0, 0), 7.0f, 1.0e-6f);
  TEST_NEAR(fields[0](0, -1, 1), 1700.0f, 1.0e-6f);
  return 0;
}
