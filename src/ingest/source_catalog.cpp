#include "gwm/ingest/source_catalog.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace gwm::ingest {

std::string to_string(SourceKind kind) {
  switch (kind) {
    case SourceKind::HRRR:
      return "HRRR";
    case SourceKind::RRFS:
      return "RRFS";
    case SourceKind::RAP:
      return "RAP";
    case SourceKind::GFS:
      return "GFS";
    case SourceKind::ERA5:
      return "ERA5";
    case SourceKind::ECMWF:
      return "ECMWF";
  }
  throw std::runtime_error("Unknown SourceKind");
}

SourceKind source_kind_from_string(const std::string& value) {
  std::string normalized = value;
  std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                 [](unsigned char ch) {
                   return static_cast<char>(std::toupper(ch));
                 });
  if (normalized == "HRRR") {
    return SourceKind::HRRR;
  }
  if (normalized == "RRFS") {
    return SourceKind::RRFS;
  }
  if (normalized == "RAP") {
    return SourceKind::RAP;
  }
  if (normalized == "GFS") {
    return SourceKind::GFS;
  }
  if (normalized == "ERA5") {
    return SourceKind::ERA5;
  }
  if (normalized == "ECMWF") {
    return SourceKind::ECMWF;
  }
  throw std::runtime_error("Unknown source string: " + value);
}

bool is_first_class_source(SourceKind kind) {
  return kind == SourceKind::HRRR || kind == SourceKind::RRFS;
}

std::unique_ptr<SourceAdapter> make_adapter(SourceKind kind) {
  switch (kind) {
    case SourceKind::HRRR:
      return std::make_unique<HrrrAdapter>();
    case SourceKind::RRFS:
      return std::make_unique<RrfsAdapter>();
    case SourceKind::RAP:
    case SourceKind::GFS:
    case SourceKind::ERA5:
    case SourceKind::ECMWF:
      return nullptr;
  }
  throw std::runtime_error("Unknown SourceKind");
}

std::vector<std::unique_ptr<SourceAdapter>> make_default_adapters() {
  std::vector<std::unique_ptr<SourceAdapter>> adapters;
  adapters.push_back(std::make_unique<HrrrAdapter>());
  adapters.push_back(std::make_unique<RrfsAdapter>());
  return adapters;
}

std::vector<ExternalToolDescriptor> recommended_external_tools() {
  return {
      {"cfrust", "C:\\Users\\drew\\cfrust", "high-level GRIB ingest", false},
      {"ecrust", "C:\\Users\\drew\\ecrust",
       "low-level GRIB/ecCodes-style access", false},
      {"wrf-rust", "C:\\Users\\drew\\wrf-rust",
       "diagnostics and verification", false},
      {"wrf-rust-plots", "C:\\Users\\drew\\wrf-rust-plots",
       "plotting and product generation", false},
      {"metrust-or-rustmet", "C:\\Users\\drew",
       "generic thermodynamic and severe-weather math", false},
      {"rusbie-or-rustbie", "C:\\Users\\drew",
       "source-data acquisition workflow", false},
  };
}

}  // namespace gwm::ingest
