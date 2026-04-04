#include "gwm/comm/collectives.hpp"

#include "gwm/comm/mpi_runtime.hpp"

#if GWM_HAVE_MPI
#include <mpi.h>
#endif

namespace gwm::comm {

void allreduce_sum_in_place(std::vector<double>& values) {
#if GWM_HAVE_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized != 0 && finalized == 0 && !values.empty()) {
    MPI_Comm comm =
        has_active_mpi_cartesian_context() ? active_mpi_cartesian_context().cart_comm
                                           : MPI_COMM_WORLD;
    MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(values.size()),
                  MPI_DOUBLE, MPI_SUM, comm);
  }
#else
  (void)values;
#endif
}

void allreduce_sum_in_place(std::vector<std::uint64_t>& values) {
#if GWM_HAVE_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized != 0 && finalized == 0 && !values.empty()) {
    MPI_Comm comm =
        has_active_mpi_cartesian_context() ? active_mpi_cartesian_context().cart_comm
                                           : MPI_COMM_WORLD;
    MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(values.size()),
                  MPI_UINT64_T, MPI_SUM, comm);
  }
#else
  (void)values;
#endif
}

}  // namespace gwm::comm
