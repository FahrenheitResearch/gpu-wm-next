#pragma once

#include <vector>

#include "gwm/core/types.hpp"

namespace gwm::surface {

class SurfaceState {
 public:
  SurfaceState(int nx, int ny, int ntile, int nsoil)
      : nx_(nx),
        ny_(ny),
        ntile_(ntile),
        nsoil_(nsoil),
        skin_temperature_(tile_cell_count(), 300.0f),
        canopy_water_(tile_cell_count(), 0.0f),
        snow_water_(tile_cell_count(), 0.0f),
        soil_temperature_(static_cast<std::size_t>(nsoil_) * tile_cell_count(),
                          295.0f),
        soil_moisture_(static_cast<std::size_t>(nsoil_) * tile_cell_count(),
                       0.25f) {}

  [[nodiscard]] int nx() const { return nx_; }
  [[nodiscard]] int ny() const { return ny_; }
  [[nodiscard]] int ntile() const { return ntile_; }
  [[nodiscard]] int nsoil() const { return nsoil_; }

  [[nodiscard]] real& skin_temperature(int tile, int i, int j) {
    return skin_temperature_[tile_index(tile, i, j)];
  }

  [[nodiscard]] real& canopy_water(int tile, int i, int j) {
    return canopy_water_[tile_index(tile, i, j)];
  }

  [[nodiscard]] real& snow_water(int tile, int i, int j) {
    return snow_water_[tile_index(tile, i, j)];
  }

  [[nodiscard]] real& soil_temperature(int layer, int tile, int i, int j) {
    return soil_temperature_[soil_index(layer, tile, i, j)];
  }

  [[nodiscard]] real& soil_moisture(int layer, int tile, int i, int j) {
    return soil_moisture_[soil_index(layer, tile, i, j)];
  }

 private:
  int nx_ = 0;
  int ny_ = 0;
  int ntile_ = 1;
  int nsoil_ = 4;

  std::vector<real> skin_temperature_;
  std::vector<real> canopy_water_;
  std::vector<real> snow_water_;
  std::vector<real> soil_temperature_;
  std::vector<real> soil_moisture_;

  [[nodiscard]] std::size_t tile_cell_count() const {
    return static_cast<std::size_t>(ntile_) * nx_ * ny_;
  }

  [[nodiscard]] std::size_t tile_index(int tile, int i, int j) const {
    return (static_cast<std::size_t>(tile) * ny_ + static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(nx_) +
           static_cast<std::size_t>(i);
  }

  [[nodiscard]] std::size_t soil_index(int layer, int tile, int i, int j) const {
    return static_cast<std::size_t>(layer) * tile_cell_count() +
           tile_index(tile, i, j);
  }
};

}  // namespace gwm::surface
