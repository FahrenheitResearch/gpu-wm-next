#include "gwm/ingest/analysis_state.hpp"
#include "gwm/ingest/boundary_cache.hpp"
#include "gwm/ingest/source_catalog.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;
  using namespace gwm::ingest;

  const auto hrrr = make_adapter(SourceKind::HRRR);
  TEST_CHECK(hrrr != nullptr);

  AnalysisStateIR analysis{};
  analysis.source = SourceKind::HRRR;
  analysis.grid.nx = 4;
  analysis.grid.ny = 3;
  analysis.grid.nz = 2;
  analysis.cycle_time_utc = "2026-04-03T00:00:00Z";
  analysis.valid_time_utc = "2026-04-03T01:00:00Z";
  analysis.forecast_offset_seconds = 3600;
  analysis.atmosphere.values["u_wind"] = std::vector<real>(analysis.grid.cell_count_3d(), 1.0f);
  analysis.atmosphere.values["v_wind"] = std::vector<real>(analysis.grid.cell_count_3d(), 2.0f);
  analysis.atmosphere.values["w_wind"] = std::vector<real>(analysis.grid.cell_count_3d(), 0.1f);
  analysis.atmosphere.values["air_temperature"] =
      std::vector<real>(analysis.grid.cell_count_3d(), 300.0f);
  analysis.atmosphere.values["water_vapor_mixing_ratio"] =
      std::vector<real>(analysis.grid.cell_count_3d(), 0.01f);
  analysis.surface.values["surface_pressure"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 95000.0f);
  analysis.surface.values["air_temperature_2m"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 299.0f);
  analysis.surface.values["specific_humidity_2m"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 0.009f);
  analysis.surface.values["u_wind_10m"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 8.0f);
  analysis.surface.values["v_wind_10m"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 1.0f);
  analysis.surface.values["skin_temperature"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 301.0f);
  analysis.static_surface.values["terrain_height"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 250.0f);
  analysis.static_surface.values["land_mask"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 1.0f);
  analysis.static_surface.values["land_use_index"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 7.0f);

  TEST_CHECK(analysis.has_required_fields(*hrrr));

  const auto canonical = analysis.to_canonical();
  TEST_CHECK(!canonical.empty());
  TEST_CHECK(canonical.atmosphere.has("u_wind"));
  TEST_CHECK(canonical.surface.has("surface_pressure"));
  TEST_CHECK(canonical.static_surface.has("terrain_height"));

  const auto round_trip = AnalysisStateIR::from_canonical(canonical);
  TEST_CHECK(round_trip.has_required_fields(*hrrr));

  BoundaryCacheIR cache{};
  cache.source = SourceKind::HRRR;
  cache.grid = analysis.grid;
  cache.cycle_time_utc = analysis.cycle_time_utc;
  cache.boundary_interval_seconds = hrrr->boundary_interval_seconds();

  BoundarySnapshotIR t0{};
  t0.forecast_offset_seconds = 0;
  t0.valid_time_utc = "2026-04-03T00:00:00Z";
  t0.atmosphere.values["u_wind"] = std::vector<real>(analysis.grid.cell_count_3d(), 1.0f);
  t0.atmosphere.values["v_wind"] = std::vector<real>(analysis.grid.cell_count_3d(), 2.0f);
  t0.atmosphere.values["w_wind"] = std::vector<real>(analysis.grid.cell_count_3d(), 0.1f);
  t0.atmosphere.values["air_temperature"] =
      std::vector<real>(analysis.grid.cell_count_3d(), 300.0f);
  t0.atmosphere.values["water_vapor_mixing_ratio"] =
      std::vector<real>(analysis.grid.cell_count_3d(), 0.01f);
  t0.surface.values["surface_pressure"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 95000.0f);
  t0.surface.values["air_temperature_2m"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 299.0f);
  t0.surface.values["specific_humidity_2m"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 0.009f);
  t0.surface.values["u_wind_10m"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 8.0f);
  t0.surface.values["v_wind_10m"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 1.0f);
  t0.surface.values["skin_temperature"] =
      std::vector<real>(analysis.grid.cell_count_2d(), 301.0f);

  auto t1 = t0;
  t1.forecast_offset_seconds = 3600;
  t1.valid_time_utc = "2026-04-03T01:00:00Z";

  cache.snapshots.push_back(t0);
  cache.snapshots.push_back(t1);

  TEST_CHECK(!cache.empty());
  TEST_CHECK(cache.has_monotonic_offsets());
  TEST_CHECK(cache.snapshots.front().has_required_fields(*hrrr));

  const auto [lower, upper] = cache.bracket(1800);
  TEST_CHECK(lower != nullptr);
  TEST_CHECK(upper != nullptr);
  TEST_CHECK(lower->forecast_offset_seconds == 0);
  TEST_CHECK(upper->forecast_offset_seconds == 3600);

  return 0;
}
