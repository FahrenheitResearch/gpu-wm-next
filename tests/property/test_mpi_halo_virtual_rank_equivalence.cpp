#include "gwm/comm/halo_exchange.hpp"
#include "gwm/comm/mpi_runtime.hpp"
#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/state/field3d.hpp"

#include "test_assert.hpp"

#if GWM_HAVE_MPI
#include <mpi.h>
#endif

namespace {

void compare_field(const gwm::state::Field3D<gwm::real>& a,
                   const gwm::state::Field3D<gwm::real>& b, double tol) {
  TEST_CHECK(a.nx() == b.nx());
  TEST_CHECK(a.ny() == b.ny());
  TEST_CHECK(a.nz() == b.nz());
  for (int k = 0; k < a.nz(); ++k) {
    for (int j = -a.halo(); j < a.ny() + a.halo(); ++j) {
      for (int i = -a.halo(); i < a.nx() + a.halo(); ++i) {
        TEST_NEAR(a(i, j, k), b(i, j, k), tol);
      }
    }
  }
}

}  // namespace

int main() {
#if GWM_HAVE_MPI
  MPI_Init(nullptr, nullptr);
  using namespace gwm;

  const auto layout =
      comm::VirtualRankLayout::build(11, 7, 3, 1, 2, 2, false, true);
  auto oracle_fields = comm::VirtualRankLayout::scatter_scalar(
      layout,
      [](int i, int j, int k) {
        return static_cast<real>(1000 * k + 100 * j + i);
      },
      "oracle");
  comm::HaloExchange::exchange_scalar(oracle_fields, layout);

  const auto context =
      comm::make_mpi_cartesian_context(MPI_COMM_WORLD, layout);
  comm::activate_mpi_cartesian_context(context);
  const auto& desc = layout[static_cast<std::size_t>(context.world_rank)];

  auto local_fields = comm::VirtualRankLayout::scatter_scalar(
      {desc},
      [](int i, int j, int k) {
        return static_cast<real>(1000 * k + 100 * j + i);
      },
      "local");

  std::vector<domain::SubdomainDescriptor> local_layout{desc};
  comm::HaloExchange::exchange_scalar(local_fields, local_layout);

  compare_field(local_fields.front(),
                oracle_fields[static_cast<std::size_t>(context.world_rank)],
                1.0e-6);

  comm::deactivate_mpi_cartesian_context();
  MPI_Finalize();
#endif
  return 0;
}
