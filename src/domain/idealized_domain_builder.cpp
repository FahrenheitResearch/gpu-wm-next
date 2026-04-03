#include "gwm/domain/idealized_domain_builder.hpp"

#include <cmath>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/core/types.hpp"

namespace gwm::domain {

namespace {

constexpr double kPi = 3.14159265358979323846;

real cosine_mountain_height(const RectilinearDomainConfig& config, int i, int j) {
  if (config.terrain_kind != TerrainProfileKind::CosineMountain ||
      config.mountain_height <= 0.0) {
    return 0.0f;
  }

  const double nxm1 = static_cast<double>(config.nx - 1);
  const double nym1 = static_cast<double>(config.ny - 1);
  const double x = (nxm1 > 0.0) ? static_cast<double>(i) / nxm1 : 0.0;
  const double y = (nym1 > 0.0) ? static_cast<double>(j) / nym1 : 0.0;

  const double dx = (x - config.mountain_center_x) / config.mountain_half_width_x;
  const double dy = (y - config.mountain_center_y) / config.mountain_half_width_y;

  if (std::abs(dx) >= 1.0 || std::abs(dy) >= 1.0) {
    return 0.0f;
  }

  const double wx = 0.5 * (1.0 + std::cos(kPi * dx));
  const double wy = 0.5 * (1.0 + std::cos(kPi * dy));
  return static_cast<real>(config.mountain_height * wx * wy);
}

}  // namespace

IdealizedDomain build_rectilinear_domain(const RectilinearDomainConfig& config) {
  gwm::require(config.nx > 0 && config.ny > 0 && config.nz > 0,
               "RectilinearDomainConfig dimensions must be positive");
  gwm::require(config.dx > 0.0 && config.dy > 0.0 && config.z_top > 0.0,
               "RectilinearDomainConfig spacing/top must be positive");
  gwm::require(config.ranks_x > 0 && config.ranks_y > 0,
               "RectilinearDomainConfig rank counts must be positive");

  IdealizedDomain domain{};
  domain.terrain_dyn.resize(static_cast<std::size_t>(config.nx) * config.ny);
  for (int j = 0; j < config.ny; ++j) {
    for (int i = 0; i < config.nx; ++i) {
      domain.terrain_dyn[static_cast<std::size_t>(j) * config.nx +
                         static_cast<std::size_t>(i)] =
          cosine_mountain_height(config, i, j);
    }
  }
  domain.metrics = GridMetrics::make_hybrid_height(
      config.nx, config.ny, config.nz, config.dx, config.dy, config.z_top,
      domain.terrain_dyn, config.terrain_taper_eta, config.periodic_x,
      config.periodic_y);
  domain.layout = comm::VirtualRankLayout::build(
      config.nx, config.ny, config.nz, config.halo, config.ranks_x,
      config.ranks_y, config.periodic_x, config.periodic_y);
  return domain;
}

}  // namespace gwm::domain
