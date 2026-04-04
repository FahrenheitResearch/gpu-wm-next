#pragma once

#include <string>
#include <vector>

#include "gwm/comm/cartesian_topology.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"

#if GWM_HAVE_MPI
#include <mpi.h>
#endif

namespace gwm::comm {

struct MpiRuntimeInfo {
  bool compiled_with_mpi = false;
  bool launcher_available = false;
  bool initialized = false;
  int world_size = 1;
  int world_rank = 0;
  std::string backend_name;
};

struct MpiCartesianContext {
  bool active = false;
  int world_size = 1;
  int world_rank = 0;
  int cart_rank = 0;
  int coord_x = 0;
  int coord_y = 0;
  int ranks_x = 1;
  int ranks_y = 1;
  bool periodic_x = true;
  bool periodic_y = true;
  CartesianNeighborRanks neighbors{};
#if GWM_HAVE_MPI
  MPI_Comm cart_comm = MPI_COMM_NULL;
#endif
};

[[nodiscard]] MpiRuntimeInfo query_mpi_runtime();

#if GWM_HAVE_MPI
[[nodiscard]] MpiCartesianContext make_mpi_cartesian_context(
    MPI_Comm parent_comm,
    const std::vector<domain::SubdomainDescriptor>& layout);
#endif

void activate_mpi_cartesian_context(MpiCartesianContext context);
void deactivate_mpi_cartesian_context();
[[nodiscard]] bool has_active_mpi_cartesian_context();
[[nodiscard]] const MpiCartesianContext& active_mpi_cartesian_context();

void release_mpi_cartesian_context(MpiCartesianContext& context);

}  // namespace gwm::comm
