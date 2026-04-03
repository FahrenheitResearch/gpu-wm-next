#pragma once

#include <string>

#include "gwm/state/field3d.hpp"

namespace gwm::state {

enum class FaceOrientation { X, Y, Z };

template <typename T>
class FaceField {
 public:
  FaceField() = default;

  FaceField(int nx, int ny, int nz, int halo, FaceOrientation orientation,
            std::string label = {})
      : orientation_(orientation) {
    int fx = nx + (orientation == FaceOrientation::X ? 1 : 0);
    int fy = ny + (orientation == FaceOrientation::Y ? 1 : 0);
    int fz = nz + (orientation == FaceOrientation::Z ? 1 : 0);
    storage_ = Field3D<T>(fx, fy, fz, halo, std::move(label));
  }

  [[nodiscard]] FaceOrientation orientation() const { return orientation_; }
  [[nodiscard]] Field3D<T>& storage() { return storage_; }
  [[nodiscard]] const Field3D<T>& storage() const { return storage_; }

 private:
  FaceOrientation orientation_ = FaceOrientation::X;
  Field3D<T> storage_;
};

}  // namespace gwm::state
