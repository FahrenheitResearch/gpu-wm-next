#include "gwm/comm/mpi_runtime.hpp"

#if GWM_HAVE_MPI
#include <mpi.h>
#endif

namespace gwm::comm {

MpiRuntimeInfo query_mpi_runtime() {
  MpiRuntimeInfo info{};
#if GWM_HAVE_MPI
  info.compiled_with_mpi = true;
  info.launcher_available = true;
  info.backend_name = "MPI";
#else
  info.compiled_with_mpi = false;
  info.launcher_available = false;
  info.backend_name = "none";
#endif
  return info;
}

}  // namespace gwm::comm
