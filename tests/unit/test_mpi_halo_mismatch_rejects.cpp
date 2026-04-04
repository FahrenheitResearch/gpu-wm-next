#include <stdexcept>

#include "gwm/comm/mpi_runtime.hpp"
#include "gwm/comm/virtual_rank_layout.hpp"

#include "test_assert.hpp"

#if GWM_HAVE_MPI
#include <mpi.h>
#endif

int main() {
#if GWM_HAVE_MPI
  MPI_Init(nullptr, nullptr);

  bool size_mismatch_threw = false;
  try {
    const auto bad_layout =
        gwm::comm::VirtualRankLayout::build(8, 8, 2, 1, 1, 1, true, true);
    auto context =
        gwm::comm::make_mpi_cartesian_context(MPI_COMM_WORLD, bad_layout);
    gwm::comm::release_mpi_cartesian_context(context);
  } catch (const std::runtime_error&) {
    size_mismatch_threw = true;
  }
  TEST_CHECK(size_mismatch_threw);

  MPI_Finalize();
#endif
  return 0;
}
