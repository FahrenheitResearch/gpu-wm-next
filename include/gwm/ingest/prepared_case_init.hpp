#pragma once

#include <string>
#include <vector>

#include "gwm/domain/grid_metrics.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/ingest/runtime_case.hpp"

namespace gwm::ingest {

struct PreparedCaseInitConfig {
  int halo = 1;
  int ranks_x = 1;
  int ranks_y = 1;
  double terrain_taper_eta = 0.25;
  bool periodic_x = false;
  bool periodic_y = false;
};

[[nodiscard]] std::vector<domain::SubdomainDescriptor> make_prepared_case_layout(
    const AnalysisStateIR& analysis, const PreparedCaseInitConfig& config);

[[nodiscard]] domain::GridMetrics make_prepared_case_metrics(
    const AnalysisStateIR& analysis, const PreparedCaseInitConfig& config);

[[nodiscard]] std::vector<dycore::DryState> make_dry_states_from_analysis(
    const AnalysisStateIR& analysis,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const std::string& label_prefix = "prepared_case");

class PreparedCaseBoundaryUpdater final : public dycore::BoundaryUpdater {
 public:
  PreparedCaseBoundaryUpdater(AnalysisStateIR analysis, BoundaryCacheIR cache,
                              PreparedCaseInitConfig config);

  void set_step_start_time(real step_start_time_seconds);

  void apply(std::vector<dycore::DryState>& states,
             const std::vector<domain::SubdomainDescriptor>& layout,
             real sim_time) override;

 private:
  AnalysisStateIR analysis_;
  BoundaryCacheIR cache_;
  PreparedCaseInitConfig config_{};
  real step_start_time_seconds_ = 0.0f;
};

}  // namespace gwm::ingest
