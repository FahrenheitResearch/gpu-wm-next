#pragma once

#include <string>
#include <vector>

#include "gwm/domain/grid_metrics.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/state/tracer_state.hpp"

namespace gwm::dycore {

class TracerBoundaryUpdater {
 public:
  virtual ~TracerBoundaryUpdater() = default;
  virtual void apply(std::vector<state::TracerState>& states,
                     const std::vector<domain::SubdomainDescriptor>& layout,
                     real sim_time) = 0;
};

class NullTracerBoundaryUpdater final : public TracerBoundaryUpdater {
 public:
  void apply(std::vector<state::TracerState>& states,
             const std::vector<domain::SubdomainDescriptor>& layout,
             real sim_time) override;
};

[[nodiscard]] std::vector<state::TracerState>
make_specific_humidity_tracers_from_global_field(
    const std::vector<DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout, int global_nx,
    int global_ny, const std::vector<real>& specific_humidity,
    const std::string& label_prefix = "specific_humidity");

void advance_passive_tracers_ssprk3(
    std::vector<state::TracerState>& tracers,
    const std::vector<DryState>& dry_states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const DryStepperConfig& config,
    TracerBoundaryUpdater& boundary_updater);

}  // namespace gwm::dycore
