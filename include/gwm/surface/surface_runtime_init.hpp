#pragma once

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "gwm/ingest/canonical_ir.hpp"
#include "gwm/surface/surface_state.hpp"
#include "gwm/surface/surface_static_properties.hpp"

namespace gwm::surface {

struct SurfaceRuntimeInitResult {
  SurfaceState state;
  SurfaceStaticProperties properties;
};

namespace detail {

[[nodiscard]] inline std::size_t surface_cell_index(int nx, int i, int j) {
  return static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

[[nodiscard]] inline const std::vector<real>& require_field(
    const gwm::ingest::FieldBundle2D& bundle, const char* field_name) {
  const auto it = bundle.values.find(field_name);
  gwm::require(it != bundle.values.end(),
               std::string("Missing canonical surface field: ") + field_name);
  return it->second;
}

inline void require_field_shape(const std::vector<real>& values, int nx, int ny,
                                const char* field_name) {
  gwm::require(static_cast<int>(values.size()) == nx * ny,
               std::string("Canonical surface field has wrong size: ") +
                   field_name);
}

[[nodiscard]] inline real cell_value(const std::vector<real>& values, int nx,
                                     int i, int j) {
  return values[surface_cell_index(nx, i, j)];
}

}  // namespace detail

[[nodiscard]] inline SurfaceState make_surface_state_from_canonical_fields(
    const gwm::ingest::FieldBundle2D& surface_fields, int nx, int ny,
    int ntile = 1, int nsoil = 4) {
  gwm::require(ntile == 1,
               "Prepared-case surface init currently supports ntile=1 only");

  SurfaceState state(nx, ny, ntile, nsoil);

  const auto& skin_temperature = detail::require_field(surface_fields,
                                                       "skin_temperature");
  detail::require_field_shape(skin_temperature, nx, ny, "skin_temperature");

  const std::vector<real>* canopy_water = nullptr;
  const std::vector<real>* snow_water = nullptr;
  const std::vector<real>* soil_temperature = nullptr;
  const std::vector<real>* soil_moisture = nullptr;

  if (surface_fields.has("canopy_water")) {
    canopy_water = &detail::require_field(surface_fields, "canopy_water");
    detail::require_field_shape(*canopy_water, nx, ny, "canopy_water");
  }
  if (surface_fields.has("snow_water")) {
    snow_water = &detail::require_field(surface_fields, "snow_water");
    detail::require_field_shape(*snow_water, nx, ny, "snow_water");
  }
  if (surface_fields.has("soil_temperature")) {
    soil_temperature = &detail::require_field(surface_fields, "soil_temperature");
    detail::require_field_shape(*soil_temperature, nx, ny, "soil_temperature");
  }
  if (surface_fields.has("soil_moisture")) {
    soil_moisture = &detail::require_field(surface_fields, "soil_moisture");
    detail::require_field_shape(*soil_moisture, nx, ny, "soil_moisture");
  }

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const real skin = detail::cell_value(skin_temperature, nx, i, j);
      state.skin_temperature(0, i, j) = skin;
      state.canopy_water(0, i, j) =
          canopy_water ? detail::cell_value(*canopy_water, nx, i, j) : 0.0f;
      state.snow_water(0, i, j) =
          snow_water ? detail::cell_value(*snow_water, nx, i, j) : 0.0f;

      const real soil_temp =
          soil_temperature ? detail::cell_value(*soil_temperature, nx, i, j)
                           : skin;
      const real soil_moist =
          soil_moisture ? detail::cell_value(*soil_moisture, nx, i, j) : 0.25f;
      for (int layer = 0; layer < nsoil; ++layer) {
        state.soil_temperature(layer, 0, i, j) = soil_temp;
        state.soil_moisture(layer, 0, i, j) = soil_moist;
      }
    }
  }

  state.validate();
  return state;
}

[[nodiscard]] inline SurfaceStaticProperties
make_surface_static_properties_from_canonical_fields(
    const gwm::ingest::FieldBundle2D& static_surface_fields, int nx, int ny,
    int ntile = 1) {
  gwm::require(ntile == 1,
               "Prepared-case surface init currently supports ntile=1 only");

  SurfaceStaticProperties properties(nx, ny, ntile);

  const auto& terrain_height =
      detail::require_field(static_surface_fields, "terrain_height");
  const auto& land_mask =
      detail::require_field(static_surface_fields, "land_mask");
  const auto& land_use_index =
      detail::require_field(static_surface_fields, "land_use_index");
  detail::require_field_shape(terrain_height, nx, ny, "terrain_height");
  detail::require_field_shape(land_mask, nx, ny, "land_mask");
  detail::require_field_shape(land_use_index, nx, ny, "land_use_index");

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const auto idx = detail::surface_cell_index(nx, i, j);
      const bool is_land = land_mask[idx] >= 0.5f;
      properties.set_tile_fraction(0, i, j, 1.0f);
      properties.set_z0m(0, i, j, is_land ? 0.10f : 0.0002f);
      properties.set_z0h(0, i, j, is_land ? 0.02f : 0.0002f);
      properties.set_kind(0, i, j,
                          is_land ? SurfaceKind::Land : SurfaceKind::Water);

      (void)terrain_height[idx];
      (void)land_use_index[idx];
    }
  }

  properties.validate();
  return properties;
}

[[nodiscard]] inline SurfaceRuntimeInitResult make_surface_runtime_from_canonical_fields(
    const gwm::ingest::FieldBundle2D& surface_fields,
    const gwm::ingest::FieldBundle2D& static_surface_fields, int nx, int ny,
    int ntile = 1, int nsoil = 4) {
  return {make_surface_state_from_canonical_fields(surface_fields, nx, ny,
                                                   ntile, nsoil),
          make_surface_static_properties_from_canonical_fields(
              static_surface_fields, nx, ny, ntile)};
}

}  // namespace gwm::surface
