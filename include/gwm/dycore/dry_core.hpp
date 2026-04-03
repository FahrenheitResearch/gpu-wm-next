#pragma once

#include <memory>
#include <string>
#include <vector>

#include "gwm/domain/grid_metrics.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"
#include "gwm/state/face_field.hpp"
#include "gwm/state/field3d.hpp"

namespace gwm::dycore {

struct DryState {
  state::Field3D<real> rho_d;
  state::Field3D<real> rho_theta_m;
  state::FaceField<real> mom_u;
  state::FaceField<real> mom_v;
  state::FaceField<real> mom_w;

  DryState() = default;
  DryState(int nx, int ny, int nz, int halo, const std::string& label_prefix);

  [[nodiscard]] DryState clone_empty_like(
      const std::string& suffix = "_tmp") const;

  void copy_all_from(const DryState& other);
  void fill_constant(real rho_value, real theta_value, real u_value,
                     real v_value, real w_value);

  [[nodiscard]] double total_dry_mass() const;
  [[nodiscard]] double total_rho_theta_m() const;
  [[nodiscard]] double total_vertical_momentum() const;
};

struct DrySlowTendencies {
  state::Field3D<real> rho_d;
  state::Field3D<real> rho_theta_m;
  state::FaceField<real> mom_u;
  state::FaceField<real> mom_v;
  state::FaceField<real> mom_w;

  DrySlowTendencies() = default;
  DrySlowTendencies(int nx, int ny, int nz, int halo,
                    const std::string& label_prefix);

  void fill_zero();
};

struct DryStepperConfig {
  real dt = 2.0f;
  int fast_substeps = 1;
  real gravity = 9.81f;
};

class BoundaryUpdater {
 public:
  virtual ~BoundaryUpdater() = default;
  virtual void apply(std::vector<DryState>& states,
                     const std::vector<domain::SubdomainDescriptor>& layout,
                     real sim_time) = 0;
};

class NullBoundaryUpdater final : public BoundaryUpdater {
 public:
  void apply(std::vector<DryState>& states,
             const std::vector<domain::SubdomainDescriptor>& layout,
             real sim_time) override;
};

class FastModeIntegrator {
 public:
  virtual ~FastModeIntegrator() = default;
  virtual void apply(std::vector<DryState>& states,
                     const std::vector<domain::SubdomainDescriptor>& layout,
                     const domain::GridMetrics& metrics,
                     const DryStepperConfig& config) = 0;
};

class LocalSplitExplicitFastMode final : public FastModeIntegrator {
 public:
  void apply(std::vector<DryState>& states,
             const std::vector<domain::SubdomainDescriptor>& layout,
             const domain::GridMetrics& metrics,
             const DryStepperConfig& config) override;
};

std::vector<DryState> make_constant_dry_state(
    const std::vector<domain::SubdomainDescriptor>& layout, real rho_value,
    real theta_value, real u_value, real v_value, real w_value,
    const std::string& label_prefix);

std::vector<DryState> make_hydrostatic_rest_state(
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, real rho_surface, real theta_ref,
    real density_scale_height, const std::string& label_prefix);

void compute_slow_tendencies(
    const std::vector<DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const DryStepperConfig& config,
    std::vector<DrySlowTendencies>& out);

void advance_dry_state_ssprk3(
    std::vector<DryState>& states,
    const std::vector<domain::SubdomainDescriptor>& layout,
    const domain::GridMetrics& metrics, const DryStepperConfig& config,
    BoundaryUpdater& boundary_updater, FastModeIntegrator& fast_modes);

}  // namespace gwm::dycore
