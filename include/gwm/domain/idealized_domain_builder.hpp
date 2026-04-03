#pragma once

#include <vector>

#include "gwm/domain/grid_metrics.hpp"
#include "gwm/domain/subdomain_descriptor.hpp"

namespace gwm::domain {

enum class TerrainProfileKind { Flat, CosineMountain };

struct RectilinearDomainConfig {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  int halo = 1;
  int ranks_x = 1;
  int ranks_y = 1;
  double dx = 1000.0;
  double dy = 1000.0;
  double z_top = 16000.0;
  TerrainProfileKind terrain_kind = TerrainProfileKind::Flat;
  double mountain_center_x = 0.5;
  double mountain_center_y = 0.5;
  double mountain_half_width_x = 0.15;
  double mountain_half_width_y = 0.15;
  double mountain_height = 0.0;
  bool periodic_x = true;
  bool periodic_y = true;
};

struct IdealizedDomain {
  GridMetrics metrics;
  std::vector<SubdomainDescriptor> layout;
  std::vector<real> terrain_dyn;

  [[nodiscard]] real terrain(int i, int j) const {
    return terrain_dyn[static_cast<std::size_t>(j) * metrics.nx +
                       static_cast<std::size_t>(i)];
  }
};

IdealizedDomain build_rectilinear_domain(const RectilinearDomainConfig& config);

}  // namespace gwm::domain
