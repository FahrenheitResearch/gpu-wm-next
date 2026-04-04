#include "gwm/comm/mpi_runtime.hpp"
#include "gwm/comm/virtual_rank_layout.hpp"

#include "test_assert.hpp"

#if GWM_HAVE_MPI
#include <mpi.h>
#endif

int main() {
#if GWM_HAVE_MPI
  MPI_Init(nullptr, nullptr);

  using namespace gwm;
  const auto layout =
      comm::VirtualRankLayout::build(8, 6, 2, 1, 2, 2, true, false);
  const auto context =
      comm::make_mpi_cartesian_context(MPI_COMM_WORLD, layout);
  const auto& desc = layout[static_cast<std::size_t>(context.world_rank)];

  TEST_CHECK(context.active);
  TEST_CHECK(context.coord_x == desc.coord_x);
  TEST_CHECK(context.coord_y == desc.coord_y);
  TEST_CHECK(context.ranks_x == desc.ranks_x);
  TEST_CHECK(context.ranks_y == desc.ranks_y);
  TEST_CHECK(context.periodic_x == desc.periodic_x);
  TEST_CHECK(context.periodic_y == desc.periodic_y);

  auto released = context;
  comm::release_mpi_cartesian_context(released);
  MPI_Finalize();
#endif
  return 0;
}
