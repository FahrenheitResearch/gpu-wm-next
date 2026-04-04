#include "gwm/ingest/runtime_case.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "gwm/ingest/json_value.hpp"
#include "gwm/ingest/source_catalog.hpp"

namespace gwm::ingest {

namespace {

SourceKind parse_source_kind(const std::string& text) {
  if (text == "HRRR") {
    return SourceKind::HRRR;
  }
  if (text == "RRFS") {
    return SourceKind::RRFS;
  }
  if (text == "RAP") {
    return SourceKind::RAP;
  }
  if (text == "GFS") {
    return SourceKind::GFS;
  }
  if (text == "ERA5") {
    return SourceKind::ERA5;
  }
  if (text == "ECMWF") {
    return SourceKind::ECMWF;
  }
  throw std::runtime_error("Unsupported source string: " + text);
}

HorizontalGridIR parse_grid(const json::Value& value) {
  const auto& object = value.as_object();
  HorizontalGridIR grid{};
  grid.nx = object.at("nx").as_int();
  grid.ny = object.at("ny").as_int();
  grid.nz = object.at("nz").as_int();
  grid.dx = object.at("dx").as_number();
  grid.dy = object.at("dy").as_number();
  if (value.has("z_top")) {
    grid.z_top = object.at("z_top").as_number();
  }
  if (value.has("ref_lat")) {
    grid.ref_lat = object.at("ref_lat").as_number();
  }
  if (value.has("ref_lon")) {
    grid.ref_lon = object.at("ref_lon").as_number();
  }
  if (value.has("truelat1")) {
    grid.truelat1 = object.at("truelat1").as_number();
  }
  if (value.has("truelat2")) {
    grid.truelat2 = object.at("truelat2").as_number();
  }
  if (value.has("stand_lon")) {
    grid.stand_lon = object.at("stand_lon").as_number();
  }
  return grid;
}

FieldBundle3D parse_field_bundle_3d(const json::Value& value) {
  FieldBundle3D bundle{};
  for (const auto& [key, entry] : value.as_object()) {
    bundle.values.emplace(key, json::to_real_vector(entry));
  }
  return bundle;
}

FieldBundle2D parse_field_bundle_2d(const json::Value& value) {
  FieldBundle2D bundle{};
  for (const auto& [key, entry] : value.as_object()) {
    bundle.values.emplace(key, json::to_real_vector(entry));
  }
  return bundle;
}

void require_bundle_sizes(const FieldBundle3D& bundle, std::size_t expected_size,
                          const char* label) {
  for (const auto& [name, values] : bundle.values) {
    gwm::require(values.size() == expected_size,
                 std::string(label) + " field has wrong size: " + name);
  }
}

void require_bundle_sizes(const FieldBundle2D& bundle, std::size_t expected_size,
                          const char* label) {
  for (const auto& [name, values] : bundle.values) {
    gwm::require(values.size() == expected_size,
                 std::string(label) + " field has wrong size: " + name);
  }
}

FieldBundle3D interpolate_bundle(const FieldBundle3D& lower,
                                 const FieldBundle3D& upper, real alpha) {
  gwm::require(lower.values.size() == upper.values.size(),
               "Boundary atmosphere bundle count mismatch");
  FieldBundle3D out{};
  for (const auto& [name, lower_values] : lower.values) {
    const auto it = upper.values.find(name);
    gwm::require(it != upper.values.end(),
                 "Upper boundary atmosphere is missing field: " + name);
    const auto& upper_values = it->second;
    gwm::require(lower_values.size() == upper_values.size(),
                 "Boundary atmosphere field size mismatch: " + name);
    auto& dest = out.values[name];
    dest.resize(lower_values.size());
    for (std::size_t idx = 0; idx < lower_values.size(); ++idx) {
      dest[idx] = (1.0f - alpha) * lower_values[idx] + alpha * upper_values[idx];
    }
  }
  return out;
}

FieldBundle2D interpolate_bundle(const FieldBundle2D& lower,
                                 const FieldBundle2D& upper, real alpha) {
  gwm::require(lower.values.size() == upper.values.size(),
               "Boundary surface bundle count mismatch");
  FieldBundle2D out{};
  for (const auto& [name, lower_values] : lower.values) {
    const auto it = upper.values.find(name);
    gwm::require(it != upper.values.end(),
                 "Upper boundary surface is missing field: " + name);
    const auto& upper_values = it->second;
    gwm::require(lower_values.size() == upper_values.size(),
                 "Boundary surface field size mismatch: " + name);
    auto& dest = out.values[name];
    dest.resize(lower_values.size());
    for (std::size_t idx = 0; idx < lower_values.size(); ++idx) {
      dest[idx] = (1.0f - alpha) * lower_values[idx] + alpha * upper_values[idx];
    }
  }
  return out;
}

}  // namespace

AnalysisStateIR load_analysis_state_json(const std::string& path) {
  const auto root = json::parse_file(path);
  gwm::require(root.at("schema_version").as_string() ==
                   "gwm-next-analysis-state/v1",
               "Unsupported analysis-state schema");

  AnalysisStateIR analysis{};
  analysis.source = parse_source_kind(root.at("source").as_string());
  analysis.grid = parse_grid(root.at("grid"));
  analysis.cycle_time_utc = root.at("cycle_time_utc").as_string();
  analysis.valid_time_utc = root.at("valid_time_utc").as_string();
  analysis.forecast_offset_seconds = root.at("forecast_offset_seconds").as_int();
  analysis.metadata = json::to_string_map(root.at("metadata"));
  analysis.atmosphere = parse_field_bundle_3d(root.at("atmosphere"));
  analysis.surface = parse_field_bundle_2d(root.at("surface"));
  analysis.static_surface = parse_field_bundle_2d(root.at("static_surface"));

  validate_runtime_analysis(analysis);
  return analysis;
}

BoundaryCacheIR load_boundary_cache_json(const std::string& path) {
  const auto root = json::parse_file(path);
  gwm::require(root.at("schema_version").as_string() ==
                   "gwm-next-boundary-cache/v1",
               "Unsupported boundary-cache schema");

  BoundaryCacheIR cache{};
  cache.source = parse_source_kind(root.at("source").as_string());
  cache.grid = parse_grid(root.at("grid"));
  cache.cycle_time_utc = root.at("cycle_time_utc").as_string();
  cache.boundary_interval_seconds = root.at("boundary_interval_seconds").as_int();
  cache.metadata = json::to_string_map(root.at("metadata"));

  for (const auto& snapshot_value : root.at("snapshots").as_array()) {
    BoundarySnapshotIR snapshot{};
    snapshot.forecast_offset_seconds =
        snapshot_value.at("forecast_offset_seconds").as_int();
    snapshot.valid_time_utc = snapshot_value.at("valid_time_utc").as_string();
    snapshot.atmosphere = parse_field_bundle_3d(snapshot_value.at("atmosphere"));
    snapshot.surface = parse_field_bundle_2d(snapshot_value.at("surface"));
    cache.snapshots.push_back(std::move(snapshot));
  }

  validate_runtime_boundary_cache(cache);
  return cache;
}

PreparedRuntimeCase load_prepared_runtime_case(
    const std::string& analysis_state_path,
    const std::string& boundary_cache_path) {
  PreparedRuntimeCase runtime_case{};
  runtime_case.analysis = load_analysis_state_json(analysis_state_path);
  runtime_case.boundary_cache = load_boundary_cache_json(boundary_cache_path);
  gwm::require(runtime_case.analysis.source == runtime_case.boundary_cache.source,
               "Analysis and boundary cache source mismatch");
  gwm::require(runtime_case.analysis.grid.nx == runtime_case.boundary_cache.grid.nx &&
                   runtime_case.analysis.grid.ny == runtime_case.boundary_cache.grid.ny &&
                   runtime_case.analysis.grid.nz == runtime_case.boundary_cache.grid.nz,
               "Analysis and boundary cache grid mismatch");
  return runtime_case;
}

BoundarySnapshotIR interpolate_boundary_snapshot(const BoundaryCacheIR& cache,
                                                 int forecast_offset_seconds) {
  validate_runtime_boundary_cache(cache);
  const auto [lower, upper] = cache.bracket(forecast_offset_seconds);
  gwm::require(lower != nullptr && upper != nullptr,
               "Boundary cache bracket returned null snapshot");

  if (lower == upper ||
      lower->forecast_offset_seconds == upper->forecast_offset_seconds) {
    return *lower;
  }

  const real alpha =
      static_cast<real>(forecast_offset_seconds - lower->forecast_offset_seconds) /
      static_cast<real>(upper->forecast_offset_seconds -
                        lower->forecast_offset_seconds);

  BoundarySnapshotIR out{};
  out.forecast_offset_seconds = forecast_offset_seconds;
  out.valid_time_utc = lower->valid_time_utc + "->" + upper->valid_time_utc;
  out.atmosphere = interpolate_bundle(lower->atmosphere, upper->atmosphere, alpha);
  out.surface = interpolate_bundle(lower->surface, upper->surface, alpha);
  return out;
}

void validate_runtime_analysis(const AnalysisStateIR& analysis) {
  gwm::require(analysis.grid.nx > 0 && analysis.grid.ny > 0 && analysis.grid.nz > 0,
               "Analysis grid dimensions must be positive");
  gwm::require(analysis.grid.dx > 0.0 && analysis.grid.dy > 0.0 &&
                   analysis.grid.z_top > 0.0,
               "Analysis grid spacing and z_top must be positive");

  const auto adapter = make_adapter(analysis.source);
  gwm::require(adapter != nullptr,
               "Analysis source is unsupported by runtime adapter");
  gwm::require(analysis.has_required_fields(*adapter),
               "Analysis state is missing required canonical fields");

  require_bundle_sizes(analysis.atmosphere, analysis.grid.cell_count_3d(),
                       "Analysis atmosphere");
  require_bundle_sizes(analysis.surface, analysis.grid.cell_count_2d(),
                       "Analysis surface");
  require_bundle_sizes(analysis.static_surface, analysis.grid.cell_count_2d(),
                       "Analysis static_surface");
}

void validate_runtime_boundary_cache(const BoundaryCacheIR& cache) {
  gwm::require(cache.grid.nx > 0 && cache.grid.ny > 0 && cache.grid.nz > 0,
               "Boundary-cache grid dimensions must be positive");
  gwm::require(cache.boundary_interval_seconds > 0,
               "Boundary-cache interval must be positive");

  const auto adapter = make_adapter(cache.source);
  gwm::require(adapter != nullptr,
               "Boundary-cache source is unsupported by runtime adapter");
  gwm::require(!cache.empty(), "Boundary cache must contain snapshots");
  gwm::require(cache.has_monotonic_offsets(),
               "Boundary cache offsets must be strictly monotonic");

  for (const auto& snapshot : cache.snapshots) {
    gwm::require(snapshot.has_required_fields(*adapter),
                 "Boundary snapshot is missing required canonical fields");
    require_bundle_sizes(snapshot.atmosphere, cache.grid.cell_count_3d(),
                         "Boundary atmosphere");
    require_bundle_sizes(snapshot.surface, cache.grid.cell_count_2d(),
                         "Boundary surface");
  }
}

}  // namespace gwm::ingest
