#pragma once

#include <string>

#include "gwm/ingest/analysis_state.hpp"
#include "gwm/ingest/boundary_cache.hpp"

namespace gwm::ingest {

struct PreparedRuntimeCase {
  AnalysisStateIR analysis;
  BoundaryCacheIR boundary_cache;
};

[[nodiscard]] AnalysisStateIR load_analysis_state_json(const std::string& path);
[[nodiscard]] BoundaryCacheIR load_boundary_cache_json(const std::string& path);
[[nodiscard]] PreparedRuntimeCase load_prepared_runtime_case(
    const std::string& analysis_state_path,
    const std::string& boundary_cache_path);

[[nodiscard]] BoundarySnapshotIR interpolate_boundary_snapshot(
    const BoundaryCacheIR& cache, int forecast_offset_seconds);

void validate_runtime_analysis(const AnalysisStateIR& analysis);
void validate_runtime_boundary_cache(const BoundaryCacheIR& cache);

}  // namespace gwm::ingest
