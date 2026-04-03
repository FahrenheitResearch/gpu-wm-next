#include "gwm/comm/mpi_runtime.hpp"

#include "test_assert.hpp"

int main() {
  const auto info = gwm::comm::query_mpi_runtime();
#if GWM_HAVE_MPI
  TEST_CHECK(info.compiled_with_mpi);
  TEST_CHECK(info.launcher_available);
  TEST_CHECK(!info.backend_name.empty());
#else
  TEST_CHECK(!info.compiled_with_mpi);
  TEST_CHECK(!info.launcher_available);
  TEST_CHECK(info.backend_name == "none");
#endif
  return 0;
}
