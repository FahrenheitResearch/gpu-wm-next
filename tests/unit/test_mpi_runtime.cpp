#include "gwm/comm/mpi_runtime.hpp"
#include "gwm/comm/virtual_rank_layout.hpp"

#include "test_assert.hpp"

int main() {
#if GWM_HAVE_MPI
  MPI_Init(nullptr, nullptr);
#endif
  const auto info = gwm::comm::query_mpi_runtime();
#if GWM_HAVE_MPI
  TEST_CHECK(info.compiled_with_mpi);
  TEST_CHECK(info.initialized);
  TEST_CHECK(!info.backend_name.empty());
  auto layout =
      gwm::comm::VirtualRankLayout::build(4, 4, 1, 1, 1, 1, true, true);
  auto context =
      gwm::comm::make_mpi_cartesian_context(MPI_COMM_WORLD, layout);
  TEST_CHECK(context.active);
  gwm::comm::release_mpi_cartesian_context(context);
  MPI_Finalize();
#else
  TEST_CHECK(!info.compiled_with_mpi);
  TEST_CHECK(!info.launcher_available);
  TEST_CHECK(info.backend_name == "none");
#endif
  return 0;
}
