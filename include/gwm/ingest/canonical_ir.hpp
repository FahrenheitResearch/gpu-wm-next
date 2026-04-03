#pragma once

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
};

struct CanonicalAnalysisIR {
  SourceKind source = SourceKind::HRRR;
  HorizontalGridIR grid;
  std::unordered_map<std::string, std::string> metadata;
  std::unordered_map<std::string, std::vector<real>> atmosphere_3d;
  std::unordered_map<std::string, std::vector<real>> surface_2d;
};

class SourceAdapter {
 public:
  virtual ~SourceAdapter() = default;
  [[nodiscard]] virtual SourceKind kind() const = 0;
  [[nodiscard]] virtual std::string name() const = 0;
  [[nodiscard]] virtual std::vector<std::string> required_fields() const = 0;
};

class HrrrAdapter final : public SourceAdapter {
 public:
  [[nodiscard]] SourceKind kind() const override { return SourceKind::HRRR; }
  [[nodiscard]] std::string name() const override { return "HRRR"; }
  [[nodiscard]] std::vector<std::string> required_fields() const override {
    return {"U", "V", "W", "T", "QVAPOR", "PSFC", "T2", "Q2", "U10", "V10"};
  }
};

class RrfsAdapter final : public SourceAdapter {
 public:
  [[nodiscard]] SourceKind kind() const override { return SourceKind::RRFS; }
  [[nodiscard]] std::string name() const override { return "RRFS"; }
  [[nodiscard]] std::vector<std::string> required_fields() const override {
    return {"UGRD", "VGRD", "DZDT", "TMP", "SPFH", "PRES", "TMP_2M", "SPFH_2M"};
  }
};

}  // namespace gwm::ingest
