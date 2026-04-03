#include "gwm/comm/halo_exchange.hpp"
#include "gwm/comm/virtual_rank_layout.hpp"

#include "test_assert.hpp"

namespace {

gwm::state::FaceField<gwm::real> make_local_face_field(
    const gwm::domain::SubdomainDescriptor& desc,
    gwm::state::FaceOrientation orientation) {
  gwm::state::FaceField<gwm::real> field(desc.nx_local(), desc.ny_local(), desc.nz,
                                         desc.halo, orientation, "face");
  field.storage().fill(-999.0f);
  for (int k = 0; k < field.storage().nz(); ++k) {
    for (int j = 0; j < field.storage().ny(); ++j) {
      for (int i = 0; i < field.storage().nx(); ++i) {
        const int gi = desc.i_begin + i;
        const int gj = desc.j_begin + j;
        field.storage()(i, j, k) = static_cast<gwm::real>(1000 * k + 100 * gj + gi);
      }
    }
  }
  return field;
}

}  // namespace

int main() {
  using namespace gwm;

  auto layout = comm::VirtualRankLayout::build(8, 8, 2, 1, 2, 2, true, true);

  std::vector<state::FaceField<real>> x_faces;
  std::vector<state::FaceField<real>> y_faces;
  std::vector<state::FaceField<real>> z_faces;
  for (const auto& desc : layout) {
    x_faces.push_back(make_local_face_field(desc, state::FaceOrientation::X));
    y_faces.push_back(make_local_face_field(desc, state::FaceOrientation::Y));
    z_faces.push_back(make_local_face_field(desc, state::FaceOrientation::Z));
  }

  comm::HaloExchange::exchange_face(x_faces, layout);
  comm::HaloExchange::exchange_face(y_faces, layout);
  comm::HaloExchange::exchange_face(z_faces, layout);

  const auto& x0 = x_faces[0].storage();
  TEST_NEAR(x0(-1, 0, 0), 7.0f, 1.0e-6f);
  TEST_NEAR(x0(5, 0, 0), 5.0f, 1.0e-6f);
  TEST_NEAR(x0(0, -1, 0), 700.0f, 1.0e-6f);
  TEST_NEAR(x0(-1, -1, 0), 707.0f, 1.0e-6f);

  const auto& y0 = y_faces[0].storage();
  TEST_NEAR(y0(-1, 0, 0), 7.0f, 1.0e-6f);
  TEST_NEAR(y0(0, -1, 0), 700.0f, 1.0e-6f);
  TEST_NEAR(y0(0, 5, 0), 500.0f, 1.0e-6f);
  TEST_NEAR(y0(-1, -1, 0), 707.0f, 1.0e-6f);

  const auto& z0 = z_faces[0].storage();
  TEST_NEAR(z0(-1, 0, 0), 7.0f, 1.0e-6f);
  TEST_NEAR(z0(0, -1, 1), 1700.0f, 1.0e-6f);
  TEST_NEAR(z0(-1, -1, 1), 1707.0f, 1.0e-6f);

  return 0;
}
