#pragma once

#include <vector>

#include "gwm/comm/cartesian_topology.hpp"

#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/state/face_field.hpp"
#include "gwm/state/field3d.hpp"

namespace gwm::comm {

struct ScalarHaloBuffers {
  std::vector<real> west_send;
  std::vector<real> east_send;
  std::vector<real> south_send;
  std::vector<real> north_send;
  std::vector<real> west_recv;
  std::vector<real> east_recv;
  std::vector<real> south_recv;
  std::vector<real> north_recv;
};

class HaloExchange {
 public:
  [[nodiscard]] static ScalarHaloBuffers allocate_scalar_buffers(
      const ScalarHaloExchangePlan& plan);

  static void pack_scalar_x_faces(const state::Field3D<real>& field,
                                  const domain::SubdomainDescriptor& desc,
                                  ScalarHaloBuffers& buffers);

  static void unpack_scalar_x_faces(state::Field3D<real>& field,
                                    const domain::SubdomainDescriptor& desc,
                                    const ScalarHaloBuffers& buffers);

  static void pack_scalar_y_faces(const state::Field3D<real>& field,
                                  const domain::SubdomainDescriptor& desc,
                                  ScalarHaloBuffers& buffers);

  static void unpack_scalar_y_faces(state::Field3D<real>& field,
                                    const domain::SubdomainDescriptor& desc,
                                    const ScalarHaloBuffers& buffers);

  static void exchange_scalar(std::vector<state::Field3D<real>>& fields,
                              const std::vector<domain::SubdomainDescriptor>& layout);

  static void exchange_face(std::vector<state::FaceField<real>>& fields,
                            const std::vector<domain::SubdomainDescriptor>& layout);

  static void exchange_face(
      const std::vector<state::FaceField<real>*>& fields,
      const std::vector<domain::SubdomainDescriptor>& layout);
};

}  // namespace gwm::comm
