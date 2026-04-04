#include "gwm/comm/mpi_runtime.hpp"

#include "gwm/core/types.hpp"

#include <utility>

#if GWM_HAVE_MPI
#include <mpi.h>
#endif

namespace gwm::comm {

namespace {

MpiCartesianContext g_active_context{};

}  // namespace

MpiRuntimeInfo query_mpi_runtime() {
  MpiRuntimeInfo info{};
#if GWM_HAVE_MPI
  info.compiled_with_mpi = true;
  info.backend_name = "MPI";

  int initialized = 0;
  MPI_Initialized(&initialized);
  info.initialized = initialized != 0;
  info.launcher_available = info.initialized;
  if (info.initialized) {
    MPI_Comm_rank(MPI_COMM_WORLD, &info.world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &info.world_size);
  }
#else
  info.compiled_with_mpi = false;
  info.launcher_available = false;
  info.initialized = false;
  info.world_size = 1;
  info.world_rank = 0;
  info.backend_name = "none";
#endif
  return info;
}

#if GWM_HAVE_MPI
MpiCartesianContext make_mpi_cartesian_context(
    MPI_Comm parent_comm, const std::vector<domain::SubdomainDescriptor>& layout) {
  gwm::require(!layout.empty(), "MPI context requires non-empty layout");

  int initialized = 0;
  MPI_Initialized(&initialized);
  gwm::require(initialized != 0,
               "MPI must be initialized before building cartesian context");

  const auto& first = layout.front();
  const int ranks_x = first.ranks_x;
  const int ranks_y = first.ranks_y;
  gwm::require(ranks_x > 0 && ranks_y > 0,
               "Cartesian context requires positive rank dimensions");
  gwm::require(static_cast<int>(layout.size()) == ranks_x * ranks_y,
               "Layout size must match ranks_x * ranks_y");

  for (const auto& desc : layout) {
    gwm::require(desc.ranks_x == ranks_x && desc.ranks_y == ranks_y,
                 "Inconsistent rank dimensions in MPI layout");
    gwm::require(desc.periodic_x == first.periodic_x &&
                     desc.periodic_y == first.periodic_y,
                 "Inconsistent periodicity in MPI layout");
  }

  int world_size = 0;
  int world_rank = 0;
  MPI_Comm_size(parent_comm, &world_size);
  MPI_Comm_rank(parent_comm, &world_rank);
  gwm::require(world_size == static_cast<int>(layout.size()),
               "MPI communicator size must match layout size");
  gwm::require(world_rank >= 0 && world_rank < world_size,
               "MPI rank must be within communicator bounds");

  const int dims[2] = {ranks_y, ranks_x};
  const int periods[2] = {first.periodic_y ? 1 : 0, first.periodic_x ? 1 : 0};
  MPI_Comm cart_comm = MPI_COMM_NULL;
  const int rc =
      MPI_Cart_create(parent_comm, 2, dims, periods, 0, &cart_comm);
  gwm::require(rc == MPI_SUCCESS && cart_comm != MPI_COMM_NULL,
               "Failed to create MPI cartesian communicator");

  int cart_rank = -1;
  MPI_Comm_rank(cart_comm, &cart_rank);
  gwm::require(cart_rank == world_rank,
               "MPI cartesian communicator reordered ranks unexpectedly");

  int coords[2] = {-1, -1};
  MPI_Cart_coords(cart_comm, cart_rank, 2, coords);
  const auto& local_desc = layout[static_cast<std::size_t>(world_rank)];
  gwm::require(local_desc.rank == world_rank,
               "Layout rank ordering must match MPI world rank");
  gwm::require(local_desc.coord_x == coords[1] && local_desc.coord_y == coords[0],
               "MPI cartesian coords must match layout coords");

  int west = MPI_PROC_NULL;
  int east = MPI_PROC_NULL;
  int south = MPI_PROC_NULL;
  int north = MPI_PROC_NULL;
  MPI_Cart_shift(cart_comm, 1, 1, &west, &east);
  MPI_Cart_shift(cart_comm, 0, 1, &south, &north);

  MpiCartesianContext context{};
  context.active = true;
  context.world_size = world_size;
  context.world_rank = world_rank;
  context.cart_rank = cart_rank;
  context.coord_x = coords[1];
  context.coord_y = coords[0];
  context.ranks_x = ranks_x;
  context.ranks_y = ranks_y;
  context.periodic_x = first.periodic_x;
  context.periodic_y = first.periodic_y;
  context.neighbors = {west == MPI_PROC_NULL ? -1 : west,
                       east == MPI_PROC_NULL ? -1 : east,
                       south == MPI_PROC_NULL ? -1 : south,
                       north == MPI_PROC_NULL ? -1 : north};
  context.cart_comm = cart_comm;

  const auto expected = find_cartesian_neighbors(layout, local_desc);
  gwm::require(expected.west == context.neighbors.west &&
                   expected.east == context.neighbors.east &&
                   expected.south == context.neighbors.south &&
                   expected.north == context.neighbors.north,
               "MPI cartesian neighbors must match layout neighbors");

  return context;
}
#endif

void activate_mpi_cartesian_context(MpiCartesianContext context) {
#if GWM_HAVE_MPI
  gwm::require(context.active, "Cannot activate an inactive MPI context");
  if (g_active_context.active) {
    release_mpi_cartesian_context(g_active_context);
  }
  g_active_context = std::move(context);
#else
  (void)context;
  gwm::require(false, "MPI context activation requested without MPI support");
#endif
}

void deactivate_mpi_cartesian_context() {
  if (g_active_context.active) {
    release_mpi_cartesian_context(g_active_context);
  }
}

bool has_active_mpi_cartesian_context() { return g_active_context.active; }

const MpiCartesianContext& active_mpi_cartesian_context() {
  gwm::require(g_active_context.active,
               "Active MPI cartesian context requested before activation");
  return g_active_context;
}

void release_mpi_cartesian_context(MpiCartesianContext& context) {
#if GWM_HAVE_MPI
  if (context.cart_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&context.cart_comm);
  }
#endif
  context = MpiCartesianContext{};
}

}  // namespace gwm::comm
