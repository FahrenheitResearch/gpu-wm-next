#include "gwm/comm/halo_exchange.hpp"
#include "gwm/comm/mpi_runtime.hpp"
#include "gwm/comm/virtual_rank_layout.hpp"

#include "test_assert.hpp"

#if GWM_HAVE_MPI
#include <mpi.h>
#endif

namespace {

gwm::state::FaceField<gwm::real> make_local_face_field(
    const gwm::domain::SubdomainDescriptor& desc,
    gwm::state::FaceOrientation orientation) {
  gwm::state::FaceField<gwm::real> field(desc.nx_local(), desc.ny_local(),
                                         desc.nz, desc.halo, orientation,
                                         "face");
  field.storage().fill(0.0f);
  for (int k = 0; k < field.storage().nz(); ++k) {
    for (int j = 0; j < field.storage().ny(); ++j) {
      for (int i = 0; i < field.storage().nx(); ++i) {
        const int gi = desc.i_begin + i;
        const int gj = desc.j_begin + j;
        field.storage()(i, j, k) =
            static_cast<gwm::real>(1000 * k + 100 * gj + gi);
      }
    }
  }
  return field;
}

void compare_face(const gwm::state::FaceField<gwm::real>& a,
                  const gwm::state::FaceField<gwm::real>& b, double tol) {
  const auto& sa = a.storage();
  const auto& sb = b.storage();
  for (int k = 0; k < sa.nz(); ++k) {
    for (int j = -sa.halo(); j < sa.ny() + sa.halo(); ++j) {
      for (int i = -sa.halo(); i < sa.nx() + sa.halo(); ++i) {
        TEST_NEAR(sa(i, j, k), sb(i, j, k), tol);
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
      comm::VirtualRankLayout::build(8, 8, 2, 1, 2, 2, true, true);

  for (const auto orientation :
       {state::FaceOrientation::X, state::FaceOrientation::Y}) {
    std::vector<state::FaceField<real>> oracle;
    for (const auto& desc : layout) {
      oracle.push_back(make_local_face_field(desc, orientation));
    }
    comm::HaloExchange::synchronize_owned_face_interfaces(oracle, layout);

    const auto context =
        comm::make_mpi_cartesian_context(MPI_COMM_WORLD, layout);
    const auto& desc = layout[static_cast<std::size_t>(context.world_rank)];
    auto local = make_local_face_field(desc, orientation);
    comm::HaloExchange::synchronize_owned_face_interfaces(local, desc,
                                                          context);

    compare_face(local, oracle[static_cast<std::size_t>(context.world_rank)],
                 1.0e-6);

    auto released = context;
    comm::release_mpi_cartesian_context(released);
  }

  MPI_Finalize();
#endif
  return 0;
}
