#pragma once

#include <utility>
#include <vector>

#include "gwm/ingest/canonical_ir.hpp"

namespace gwm::ingest {

struct BoundarySnapshotIR {
  int forecast_offset_seconds = 0;
  std::string valid_time_utc;
  FieldBundle3D atmosphere;
  FieldBundle2D surface;

  [[nodiscard]] bool has_required_fields(const SourceAdapter& adapter) const {
    for (const auto& binding : adapter.atmosphere_field_bindings()) {
      if (!atmosphere.has(binding.canonical_name)) {
        return false;
      }
    }
    for (const auto& binding : adapter.surface_field_bindings()) {
      if (!surface.has(binding.canonical_name)) {
        return false;
      }
    }
    return true;
  }
};

struct BoundaryCacheIR {
  SourceKind source = SourceKind::HRRR;
  HorizontalGridIR grid;
  std::string cycle_time_utc;
  int boundary_interval_seconds = 3600;
  std::unordered_map<std::string, std::string> metadata;
  std::vector<BoundarySnapshotIR> snapshots;

  [[nodiscard]] bool empty() const { return snapshots.empty(); }

  [[nodiscard]] bool has_monotonic_offsets() const {
    for (std::size_t idx = 1; idx < snapshots.size(); ++idx) {
      if (snapshots[idx - 1].forecast_offset_seconds >=
          snapshots[idx].forecast_offset_seconds) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] std::pair<const BoundarySnapshotIR*, const BoundarySnapshotIR*>
  bracket(int forecast_offset_seconds) const {
    if (snapshots.empty()) {
      return {nullptr, nullptr};
    }

    const BoundarySnapshotIR* lower = nullptr;
    const BoundarySnapshotIR* upper = nullptr;
    for (const auto& snapshot : snapshots) {
      if (snapshot.forecast_offset_seconds <= forecast_offset_seconds) {
        lower = &snapshot;
      }
      if (snapshot.forecast_offset_seconds >= forecast_offset_seconds) {
        upper = &snapshot;
        break;
      }
    }

    if (lower == nullptr) {
      lower = &snapshots.front();
    }
    if (upper == nullptr) {
      upper = &snapshots.back();
    }
    return {lower, upper};
  }
};

}  // namespace gwm::ingest
