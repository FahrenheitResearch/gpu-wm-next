#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "gwm/core/types.hpp"

namespace gwm::surface {

class SurfaceState {
 public:
  SurfaceState(int nx, int ny, int ntile, int nsoil)
      : nx_(require_dimension(nx, "nx")),
        ny_(require_dimension(ny, "ny")),
        ntile_(require_dimension(ntile, "ntile")),
        nsoil_(require_dimension(nsoil, "nsoil")),
        skin_temperature_(tile_cell_count(), 300.0f),
        canopy_water_(tile_cell_count(), 0.0f),
        snow_water_(tile_cell_count(), 0.0f),
        soil_temperature_(static_cast<std::size_t>(nsoil_) * tile_cell_count(),
                          295.0f),
        soil_moisture_(static_cast<std::size_t>(nsoil_) * tile_cell_count(),
                       0.25f) {
    validate();
  }

  [[nodiscard]] int nx() const { return nx_; }
  [[nodiscard]] int ny() const { return ny_; }
  [[nodiscard]] int ntile() const { return ntile_; }
  [[nodiscard]] int nsoil() const { return nsoil_; }

  [[nodiscard]] std::size_t cell_count() const {
    return static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_);
  }

  [[nodiscard]] std::size_t tile_cell_count() const {
    return static_cast<std::size_t>(ntile_) * cell_count();
  }

  [[nodiscard]] std::size_t soil_storage_count() const {
    return static_cast<std::size_t>(nsoil_) * tile_cell_count();
  }

  [[nodiscard]] std::size_t tile_linear_index(int tile, int i, int j) const {
    return tile_index(tile, i, j);
  }

  [[nodiscard]] std::size_t soil_linear_index(int layer, int tile, int i,
                                              int j) const {
    return soil_index(layer, tile, i, j);
  }

  [[nodiscard]] real& skin_temperature(int tile, int i, int j) {
    return skin_temperature_[tile_index(tile, i, j)];
  }
  [[nodiscard]] const real& skin_temperature(int tile, int i, int j) const {
    return skin_temperature_[tile_index(tile, i, j)];
  }

  [[nodiscard]] real& canopy_water(int tile, int i, int j) {
    return canopy_water_[tile_index(tile, i, j)];
  }
  [[nodiscard]] const real& canopy_water(int tile, int i, int j) const {
    return canopy_water_[tile_index(tile, i, j)];
  }

  [[nodiscard]] real& snow_water(int tile, int i, int j) {
    return snow_water_[tile_index(tile, i, j)];
  }
  [[nodiscard]] const real& snow_water(int tile, int i, int j) const {
    return snow_water_[tile_index(tile, i, j)];
  }

  [[nodiscard]] real& soil_temperature(int layer, int tile, int i, int j) {
    return soil_temperature_[soil_index(layer, tile, i, j)];
  }
  [[nodiscard]] const real& soil_temperature(int layer, int tile, int i,
                                             int j) const {
    return soil_temperature_[soil_index(layer, tile, i, j)];
  }

  [[nodiscard]] real& soil_moisture(int layer, int tile, int i, int j) {
    return soil_moisture_[soil_index(layer, tile, i, j)];
  }
  [[nodiscard]] const real& soil_moisture(int layer, int tile, int i,
                                          int j) const {
    return soil_moisture_[soil_index(layer, tile, i, j)];
  }

  void validate() const {
    gwm::require(nx_ > 0 && ny_ > 0 && ntile_ > 0 && nsoil_ > 0,
                 "SurfaceState dimensions must remain positive");
    const auto expected_tile_cells = tile_cell_count();
    const auto expected_soil_cells = soil_storage_count();
    gwm::require(skin_temperature_.size() == expected_tile_cells,
                 "SurfaceState skin_temperature storage size mismatch");
    gwm::require(canopy_water_.size() == expected_tile_cells,
                 "SurfaceState canopy_water storage size mismatch");
    gwm::require(snow_water_.size() == expected_tile_cells,
                 "SurfaceState snow_water storage size mismatch");
    gwm::require(soil_temperature_.size() == expected_soil_cells,
                 "SurfaceState soil_temperature storage size mismatch");
    gwm::require(soil_moisture_.size() == expected_soil_cells,
                 "SurfaceState soil_moisture storage size mismatch");
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

  [[nodiscard]] static int require_dimension(int value, const char* label) {
    gwm::require(value > 0, std::string("SurfaceState requires positive ") + label);
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

  [[nodiscard]] std::size_t soil_index(int layer, int tile, int i, int j) const {
    validate_soil_layer(layer);
    return static_cast<std::size_t>(layer) * tile_cell_count() +
           tile_index(tile, i, j);
  }

  void validate_tile(int tile) const {
    gwm::require(tile >= 0 && tile < ntile_, "SurfaceState tile index out of range");
  }

  void validate_soil_layer(int layer) const {
    gwm::require(layer >= 0 && layer < nsoil_,
                 "SurfaceState soil layer index out of range");
  }

  void validate_cell(int i, int j) const {
    gwm::require(i >= 0 && i < nx_ && j >= 0 && j < ny_,
                 "SurfaceState cell index out of range");
  }
};

}  // namespace gwm::surface
