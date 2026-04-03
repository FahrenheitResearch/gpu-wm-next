#pragma once

#include <algorithm>
#include <string>

#include "gwm/core/cuda_utils.hpp"
#include "gwm/core/types.hpp"

namespace gwm::state {

template <typename T>
class Field3D {
 public:
  Field3D() = default;

  Field3D(int nx, int ny, int nz, int halo, std::string label = {})
      : nx_(nx), ny_(ny), nz_(nz), halo_(halo), label_(std::move(label)) {
    allocate();
  }

  Field3D(const Field3D&) = delete;
  Field3D& operator=(const Field3D&) = delete;

  Field3D(Field3D&& other) noexcept { move_from(std::move(other)); }

  Field3D& operator=(Field3D&& other) noexcept {
    if (this != &other) {
      reset();
      move_from(std::move(other));
    }
    return *this;
  }

  ~Field3D() { reset(); }

  [[nodiscard]] int nx() const { return nx_; }
  [[nodiscard]] int ny() const { return ny_; }
  [[nodiscard]] int nz() const { return nz_; }
  [[nodiscard]] int halo() const { return halo_; }
  [[nodiscard]] int pitch_x() const { return nx_ + 2 * halo_; }
  [[nodiscard]] int pitch_y() const { return ny_ + 2 * halo_; }
  [[nodiscard]] std::size_t total_size() const {
    return static_cast<std::size_t>(pitch_x()) * pitch_y() * nz_;
  }

  [[nodiscard]] T* data() { return data_; }
  [[nodiscard]] const T* data() const { return data_; }

  [[nodiscard]] std::size_t storage_index(int i, int j, int k) const {
    return (static_cast<std::size_t>(k) * pitch_y() +
            static_cast<std::size_t>(j + halo_)) *
               static_cast<std::size_t>(pitch_x()) +
           static_cast<std::size_t>(i + halo_);
  }

  [[nodiscard]] T& operator()(int i, int j, int k) {
    return data_[storage_index(i, j, k)];
  }

  [[nodiscard]] const T& operator()(int i, int j, int k) const {
    return data_[storage_index(i, j, k)];
  }

  void fill(T value) {
    for (std::size_t idx = 0; idx < total_size(); ++idx) {
      data_[idx] = value;
    }
  }

  void copy_all_from(const Field3D& other) {
    gwm::require(nx_ == other.nx_ && ny_ == other.ny_ && nz_ == other.nz_ &&
                     halo_ == other.halo_,
                 "Field3D shape mismatch in copy_all_from");
    std::copy(other.data_, other.data_ + other.total_size(), data_);
  }

  [[nodiscard]] Field3D clone_empty_like(
      const std::string& label_suffix = "_tmp") const {
    return Field3D(nx_, ny_, nz_, halo_, label_ + label_suffix);
  }

  [[nodiscard]] double owned_sum() const {
    double sum = 0.0;
    for (int k = 0; k < nz_; ++k) {
      for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
          sum += static_cast<double>((*this)(i, j, k));
        }
      }
    }
    return sum;
  }

  void sync() const { GWM_CUDA_CHECK(cudaDeviceSynchronize()); }

 private:
  int nx_ = 0;
  int ny_ = 0;
  int nz_ = 0;
  int halo_ = 0;
  T* data_ = nullptr;
  std::string label_;

  void allocate() {
    if (total_size() == 0) {
      return;
    }
    GWM_CUDA_CHECK(cudaMallocManaged(&data_, total_size() * sizeof(T)));
    fill(T{});
  }

  void reset() noexcept {
    if (data_ != nullptr) {
      (void)cudaFree(data_);
      data_ = nullptr;
    }
  }

  void move_from(Field3D&& other) {
    nx_ = other.nx_;
    ny_ = other.ny_;
    nz_ = other.nz_;
    halo_ = other.halo_;
    data_ = other.data_;
    label_ = std::move(other.label_);
    other.nx_ = 0;
    other.ny_ = 0;
    other.nz_ = 0;
    other.halo_ = 0;
    other.data_ = nullptr;
  }
};

}  // namespace gwm::state
