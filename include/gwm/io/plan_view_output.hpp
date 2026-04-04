#pragma once

#include <string>
#include <vector>

#include "gwm/core/types.hpp"
#include "gwm/domain/grid_metrics.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/dycore/dry_core.hpp"

namespace gwm::io {

struct PlanViewField {
  std::string name;
  std::string units;
  std::string location;
  int nx = 0;
  int ny = 0;
  std::vector<double> values;
};

struct PlanViewBundle {
  std::string schema_version = "gwm-next-plan-view/v1";
  std::string case_kind;
  int steps = 0;
  real dt = 0.0f;
  int slice_k = 0;
  int nx = 0;
  int ny = 0;
  double dx = 0.0;
  double dy = 0.0;
  double slice_mean_height_m = 0.0;
  std::vector<PlanViewField> fields;
};

[[nodiscard]] PlanViewBundle extract_dry_plan_view(
    const std::vector<dycore::DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const std::string& case_kind, int steps,
    real dt, int slice_k);

[[nodiscard]] std::string plan_view_bundle_to_json(
    const PlanViewBundle& bundle);

void write_plan_view_bundle_json(const PlanViewBundle& bundle,
                                 const std::string& path);

}  // namespace gwm::io
