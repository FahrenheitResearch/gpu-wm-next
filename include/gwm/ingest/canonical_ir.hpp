#pragma once

#include <algorithm>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include "gwm/core/types.hpp"

namespace gwm::ingest {

enum class SourceKind { HRRR, RRFS, RAP, GFS, ERA5, ECMWF };

struct HorizontalGridIR {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  double dx = 0.0;
  double dy = 0.0;
  double ref_lat = 0.0;
  double ref_lon = 0.0;
  double truelat1 = 0.0;
  double truelat2 = 0.0;
  double stand_lon = 0.0;

  [[nodiscard]] std::size_t cell_count_2d() const {
    return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
  }

  [[nodiscard]] std::size_t cell_count_3d() const {
    return cell_count_2d() * static_cast<std::size_t>(nz);
  }
};

struct FieldBundle3D {
  std::unordered_map<std::string, std::vector<real>> values;

  [[nodiscard]] bool has(const std::string& name) const {
    return values.find(name) != values.end();
  }

  [[nodiscard]] std::size_t field_count() const { return values.size(); }
};

struct FieldBundle2D {
  std::unordered_map<std::string, std::vector<real>> values;

  [[nodiscard]] bool has(const std::string& name) const {
    return values.find(name) != values.end();
  }

  [[nodiscard]] std::size_t field_count() const { return values.size(); }
};

struct SourceFieldBinding {
  std::string source_name;
  std::string canonical_name;
};

struct CanonicalAnalysisIR {
  SourceKind source = SourceKind::HRRR;
  HorizontalGridIR grid;
  std::string cycle_time_utc;
  std::string valid_time_utc;
  int forecast_offset_seconds = 0;
  std::unordered_map<std::string, std::string> metadata;
  FieldBundle3D atmosphere;
  FieldBundle2D surface;
  FieldBundle2D static_surface;

  [[nodiscard]] bool empty() const {
    return atmosphere.values.empty() && surface.values.empty() &&
           static_surface.values.empty();
  }
};

class SourceAdapter {
 public:
  virtual ~SourceAdapter() = default;
  [[nodiscard]] virtual SourceKind kind() const = 0;
  [[nodiscard]] virtual std::string name() const = 0;
  [[nodiscard]] virtual std::vector<SourceFieldBinding> atmosphere_field_bindings()
      const = 0;
  [[nodiscard]] virtual std::vector<SourceFieldBinding> surface_field_bindings()
      const = 0;
  [[nodiscard]] virtual std::vector<SourceFieldBinding> static_field_bindings()
      const = 0;
  [[nodiscard]] virtual int boundary_interval_seconds() const = 0;

  [[nodiscard]] virtual std::vector<std::string> required_fields() const {
    std::vector<std::string> fields;
    const auto append_unique = [&fields](const auto& bindings) {
      for (const auto& binding : bindings) {
        const auto found = std::find(fields.begin(), fields.end(),
                                     binding.source_name);
        if (found == fields.end()) {
          fields.push_back(binding.source_name);
        }
      }
    };
    append_unique(atmosphere_field_bindings());
    append_unique(surface_field_bindings());
    append_unique(static_field_bindings());
    return fields;
  }
};

class HrrrAdapter final : public SourceAdapter {
 public:
  [[nodiscard]] SourceKind kind() const override { return SourceKind::HRRR; }
  [[nodiscard]] std::string name() const override { return "HRRR"; }
  [[nodiscard]] std::vector<SourceFieldBinding> atmosphere_field_bindings()
      const override {
    return {{"U", "u_wind"},
            {"V", "v_wind"},
            {"W", "w_wind"},
            {"T", "air_temperature"},
            {"QVAPOR", "water_vapor_mixing_ratio"}};
  }

  [[nodiscard]] std::vector<SourceFieldBinding> surface_field_bindings()
      const override {
    return {{"PSFC", "surface_pressure"},
            {"T2", "air_temperature_2m"},
            {"Q2", "specific_humidity_2m"},
            {"U10", "u_wind_10m"},
            {"V10", "v_wind_10m"},
            {"TSK", "skin_temperature"}};
  }

  [[nodiscard]] std::vector<SourceFieldBinding> static_field_bindings()
      const override {
    return {{"HGT", "terrain_height"},
            {"LANDMASK", "land_mask"},
            {"LU_INDEX", "land_use_index"}};
  }

  [[nodiscard]] int boundary_interval_seconds() const override {
    return 3600;
  }
};

class RrfsAdapter final : public SourceAdapter {
 public:
  [[nodiscard]] SourceKind kind() const override { return SourceKind::RRFS; }
  [[nodiscard]] std::string name() const override { return "RRFS"; }
  [[nodiscard]] std::vector<SourceFieldBinding> atmosphere_field_bindings()
      const override {
    return {{"UGRD", "u_wind"},
            {"VGRD", "v_wind"},
            {"DZDT", "w_wind"},
            {"TMP", "air_temperature"},
            {"SPFH", "specific_humidity"}};
  }

  [[nodiscard]] std::vector<SourceFieldBinding> surface_field_bindings()
      const override {
    return {{"PRES", "surface_pressure"},
            {"TMP_2M", "air_temperature_2m"},
            {"SPFH_2M", "specific_humidity_2m"},
            {"UGRD_10M", "u_wind_10m"},
            {"VGRD_10M", "v_wind_10m"},
            {"TMP_SFC", "skin_temperature"}};
  }

  [[nodiscard]] std::vector<SourceFieldBinding> static_field_bindings()
      const override {
    return {{"HGT_SFC", "terrain_height"},
            {"LANDMASK", "land_mask"},
            {"LANDUSE", "land_use_index"}};
  }

  [[nodiscard]] int boundary_interval_seconds() const override {
    return 3600;
  }
};

}  // namespace gwm::ingest
