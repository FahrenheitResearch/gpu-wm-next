#pragma once

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "gwm/core/types.hpp"

namespace gwm::surface {

enum class SurfaceKind { Land, Water, Snow, Ice };

class SurfaceStaticProperties {
 public:
  SurfaceStaticProperties(int nx, int ny, int ntile)
      : nx_(require_dimension(nx, "nx")),
        ny_(require_dimension(ny, "ny")),
        ntile_(require_dimension(ntile, "ntile")),
        tile_fraction_(tile_cell_count(), 1.0f / static_cast<real>(ntile_)),
        z0m_(tile_cell_count(), 0.1f),
        z0h_(tile_cell_count(), 0.02f),
        kind_(tile_cell_count(), SurfaceKind::Land) {
    validate();
  }

  [[nodiscard]] int nx() const { return nx_; }
  [[nodiscard]] int ny() const { return ny_; }
  [[nodiscard]] int ntile() const { return ntile_; }
  [[nodiscard]] bool is_single_tile() const { return ntile_ == 1; }

  [[nodiscard]] std::size_t cell_count() const {
    return static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_);
  }

  [[nodiscard]] std::size_t tile_cell_count() const {
    return static_cast<std::size_t>(ntile_) * cell_count();
  }

  [[nodiscard]] std::size_t tile_linear_index(int tile, int i, int j) const {
    return tile_index(tile, i, j);
  }

  [[nodiscard]] real& tile_fraction(int tile, int i, int j) {
    return tile_fraction_.at(tile_index(tile, i, j));
  }
  [[nodiscard]] real tile_fraction(int tile, int i, int j) const {
    return tile_fraction_.at(tile_index(tile, i, j));
  }
  void set_tile_fraction(int tile, int i, int j, real value) {
    gwm::require(value >= 0.0f,
                 "SurfaceStaticProperties tile fractions must be nonnegative");
    tile_fraction(tile, i, j) = value;
  }

  [[nodiscard]] real& roughness_momentum(int tile, int i, int j) {
    return z0m_.at(tile_index(tile, i, j));
  }
  [[nodiscard]] real roughness_momentum(int tile, int i, int j) const {
    return z0m_.at(tile_index(tile, i, j));
  }
  [[nodiscard]] real z0m(int tile, int i, int j) const {
    return roughness_momentum(tile, i, j);
  }
  void set_z0m(int tile, int i, int j, real value) {
    gwm::require(value > 0.0f,
                 "SurfaceStaticProperties momentum roughness must be positive");
    roughness_momentum(tile, i, j) = value;
  }

  [[nodiscard]] real& roughness_heat(int tile, int i, int j) {
    return z0h_.at(tile_index(tile, i, j));
  }
  [[nodiscard]] real roughness_heat(int tile, int i, int j) const {
    return z0h_.at(tile_index(tile, i, j));
  }
  [[nodiscard]] real z0h(int tile, int i, int j) const {
    return roughness_heat(tile, i, j);
  }
  void set_z0h(int tile, int i, int j, real value) {
    gwm::require(value > 0.0f,
                 "SurfaceStaticProperties heat roughness must be positive");
    roughness_heat(tile, i, j) = value;
  }

  [[nodiscard]] SurfaceKind& kind(int tile, int i, int j) {
    return kind_.at(tile_index(tile, i, j));
  }
  [[nodiscard]] SurfaceKind kind(int tile, int i, int j) const {
    return kind_.at(tile_index(tile, i, j));
  }
  void set_kind(int tile, int i, int j, SurfaceKind value) {
    kind(tile, i, j) = value;
  }

  [[nodiscard]] real tile_fraction_sum(int i, int j) const {
    validate_cell(i, j);
    real sum = 0.0f;
    for (int tile = 0; tile < ntile_; ++tile) {
      sum += tile_fraction(tile, i, j);
    }
    return sum;
  }

  void validate_cell_tile_fractions(int i, int j,
                                    real tolerance = 1.0e-5f) const {
    const real sum = tile_fraction_sum(i, j);
    gwm::require(sum > 0.0f,
                 "SurfaceStaticProperties tile fractions must sum positive");
    gwm::require(std::abs(sum - 1.0f) <= tolerance,
                 "SurfaceStaticProperties tile fractions must sum to one");
  }

  void validate(real tolerance = 1.0e-5f) const {
    const auto expected_tile_cells = tile_cell_count();
    gwm::require(tile_fraction_.size() == expected_tile_cells,
                 "SurfaceStaticProperties tile_fraction storage size mismatch");
    gwm::require(z0m_.size() == expected_tile_cells,
                 "SurfaceStaticProperties z0m storage size mismatch");
    gwm::require(z0h_.size() == expected_tile_cells,
                 "SurfaceStaticProperties z0h storage size mismatch");
    gwm::require(kind_.size() == expected_tile_cells,
                 "SurfaceStaticProperties kind storage size mismatch");

    for (int j = 0; j < ny_; ++j) {
      for (int i = 0; i < nx_; ++i) {
        validate_cell_tile_fractions(i, j, tolerance);
        for (int tile = 0; tile < ntile_; ++tile) {
          gwm::require(roughness_momentum(tile, i, j) > 0.0f,
                       "SurfaceStaticProperties momentum roughness must be "
                       "positive");
          gwm::require(
              roughness_heat(tile, i, j) > 0.0f,
              "SurfaceStaticProperties heat roughness must be positive");
        }
      }
    }
  }

 private:
  int nx_ = 0;
  int ny_ = 0;
  int ntile_ = 1;

  std::vector<real> tile_fraction_;
  std::vector<real> z0m_;
  std::vector<real> z0h_;
  std::vector<SurfaceKind> kind_;

  [[nodiscard]] static int require_dimension(int value, const char* label) {
    gwm::require(
        value > 0,
        std::string("SurfaceStaticProperties requires positive ") + label);
    return value;
  }

  [[nodiscard]] std::size_t tile_index(int tile, int i, int j) const {
    validate_tile(tile);
    validate_cell(i, j);
    return (static_cast<std::size_t>(tile) * static_cast<std::size_t>(ny_) +
            static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(nx_) +
           static_cast<std::size_t>(i);
  }

  void validate_tile(int tile) const {
    gwm::require(tile >= 0 && tile < ntile_,
                 "SurfaceStaticProperties tile index out of range");
  }

  void validate_cell(int i, int j) const {
    gwm::require(i >= 0 && i < nx_ && j >= 0 && j < ny_,
                 "SurfaceStaticProperties cell index out of range");
  }
};

}  // namespace gwm::surface
