#pragma once

#include <string>

#include "gwm/ingest/canonical_ir.hpp"

namespace gwm::ingest {

struct AnalysisStateIR {
  SourceKind source = SourceKind::HRRR;
  HorizontalGridIR grid;
  std::string cycle_time_utc;
  std::string valid_time_utc;
  int forecast_offset_seconds = 0;
  std::unordered_map<std::string, std::string> metadata;
  FieldBundle3D atmosphere;
  FieldBundle2D surface;
  FieldBundle2D static_surface;

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
    for (const auto& binding : adapter.static_field_bindings()) {
      if (!static_surface.has(binding.canonical_name)) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] CanonicalAnalysisIR to_canonical() const {
    CanonicalAnalysisIR ir{};
    ir.source = source;
    ir.grid = grid;
    ir.cycle_time_utc = cycle_time_utc;
    ir.valid_time_utc = valid_time_utc;
    ir.forecast_offset_seconds = forecast_offset_seconds;
    ir.metadata = metadata;
    ir.atmosphere = atmosphere;
    ir.surface = surface;
    ir.static_surface = static_surface;
    return ir;
  }

  [[nodiscard]] static AnalysisStateIR from_canonical(
      const CanonicalAnalysisIR& ir) {
    AnalysisStateIR state{};
    state.source = ir.source;
    state.grid = ir.grid;
    state.cycle_time_utc = ir.cycle_time_utc;
    state.valid_time_utc = ir.valid_time_utc;
    state.forecast_offset_seconds = ir.forecast_offset_seconds;
    state.metadata = ir.metadata;
    state.atmosphere = ir.atmosphere;
    state.surface = ir.surface;
    state.static_surface = ir.static_surface;
    return state;
  }
};

}  // namespace gwm::ingest
