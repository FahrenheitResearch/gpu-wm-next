#pragma once

#include <vector>

namespace gwm::domain {

struct GridMetrics {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  double dx = 0.0;
  double dy = 0.0;
  double z_top = 0.0;
  double dz_nominal = 0.0;
  double terrain_taper_eta = 0.25;
  std::vector<double> eta_interfaces;
  std::vector<double> eta_centers;
  std::vector<double> z_interfaces;
  std::vector<double> z_centers;

  static GridMetrics make_hybrid_height(int nx, int ny, int nz, double dx,
                                        double dy, double z_top) {
    GridMetrics metrics{};
    metrics.nx = nx;
    metrics.ny = ny;
    metrics.nz = nz;
    metrics.dx = dx;
    metrics.dy = dy;
    metrics.z_top = z_top;
    metrics.dz_nominal = z_top / static_cast<double>(nz);
    metrics.eta_interfaces.resize(static_cast<std::size_t>(nz + 1));
    metrics.eta_centers.resize(static_cast<std::size_t>(nz));
    metrics.z_interfaces.resize(static_cast<std::size_t>(nz + 1));
    metrics.z_centers.resize(static_cast<std::size_t>(nz));
    for (int k = 0; k <= nz; ++k) {
      metrics.eta_interfaces[static_cast<std::size_t>(k)] =
          static_cast<double>(k) / static_cast<double>(nz);
      metrics.z_interfaces[static_cast<std::size_t>(k)] =
          metrics.eta_interfaces[static_cast<std::size_t>(k)] * z_top;
      if (k < nz) {
        metrics.eta_centers[static_cast<std::size_t>(k)] =
            (static_cast<double>(k) + 0.5) / static_cast<double>(nz);
        metrics.z_centers[static_cast<std::size_t>(k)] =
            metrics.eta_centers[static_cast<std::size_t>(k)] * z_top;
      }
    }
    return metrics;
  }
};

}  // namespace gwm::domain
