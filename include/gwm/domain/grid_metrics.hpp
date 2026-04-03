#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "gwm/core/types.hpp"
#include "gwm/state/field3d.hpp"

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
  bool periodic_x = true;
  bool periodic_y = true;
  std::vector<double> eta_interfaces;
  std::vector<double> eta_centers;
  std::vector<double> terrain_weights_interfaces;
  std::vector<double> terrain_weights_centers;
  std::vector<double> z_interfaces;
  std::vector<double> z_centers;
  state::Field3D<real> terrain;
  state::Field3D<real> terrain_slope_x;
  state::Field3D<real> terrain_slope_y;
  state::Field3D<real> z_interfaces_metric;
  state::Field3D<real> z_centers_metric;
  state::Field3D<real> inv_dz_cell_metric;
  state::Field3D<real> inv_dz_face_metric;

  GridMetrics() = default;
  GridMetrics(const GridMetrics&) = delete;
  GridMetrics& operator=(const GridMetrics&) = delete;
  GridMetrics(GridMetrics&&) noexcept = default;
  GridMetrics& operator=(GridMetrics&&) noexcept = default;

  [[nodiscard]] real terrain_height(int i, int j) const {
    return terrain(i, j, 0);
  }

  [[nodiscard]] real terrain_slope_x_value(int i, int j) const {
    return terrain_slope_x(i, j, 0);
  }

  [[nodiscard]] real terrain_slope_y_value(int i, int j) const {
    return terrain_slope_y(i, j, 0);
  }

  [[nodiscard]] real z_interface(int i, int j, int k) const {
    return z_interfaces_metric(i, j, k);
  }

  [[nodiscard]] real z_center(int i, int j, int k) const {
    return z_centers_metric(i, j, k);
  }

  [[nodiscard]] real inv_dz_cell(int i, int j, int k) const {
    return inv_dz_cell_metric(i, j, k);
  }

  [[nodiscard]] real inv_dz_face(int i, int j, int k) const {
    return inv_dz_face_metric(i, j, k);
  }

  [[nodiscard]] double terrain_weight_interface(int k) const {
    return terrain_weights_interfaces[static_cast<std::size_t>(k)];
  }

  [[nodiscard]] double terrain_weight_center(int k) const {
    return terrain_weights_centers[static_cast<std::size_t>(k)];
  }

  [[nodiscard]] double flat_z_interface(int k) const {
    return z_interfaces[static_cast<std::size_t>(k)];
  }

  [[nodiscard]] double flat_z_center(int k) const {
    return z_centers[static_cast<std::size_t>(k)];
  }

  static GridMetrics make_hybrid_height(int nx, int ny, int nz, double dx,
                                        double dy, double z_top,
                                        const std::vector<real>& terrain_dyn,
                                        double terrain_taper_eta = 0.25,
                                        bool periodic_x = true,
                                        bool periodic_y = true) {
    GridMetrics metrics{};
    gwm::require(static_cast<int>(terrain_dyn.size()) == nx * ny,
                 "terrain_dyn shape mismatch in GridMetrics::make_hybrid_height");
    metrics.nx = nx;
    metrics.ny = ny;
    metrics.nz = nz;
    metrics.dx = dx;
    metrics.dy = dy;
    metrics.z_top = z_top;
    metrics.dz_nominal = z_top / static_cast<double>(nz);
    metrics.terrain_taper_eta = terrain_taper_eta;
    metrics.periodic_x = periodic_x;
    metrics.periodic_y = periodic_y;
    metrics.eta_interfaces.resize(static_cast<std::size_t>(nz + 1));
    metrics.eta_centers.resize(static_cast<std::size_t>(nz));
    metrics.terrain_weights_interfaces.resize(static_cast<std::size_t>(nz + 1));
    metrics.terrain_weights_centers.resize(static_cast<std::size_t>(nz));
    metrics.z_interfaces.resize(static_cast<std::size_t>(nz + 1));
    metrics.z_centers.resize(static_cast<std::size_t>(nz));
    metrics.terrain = state::Field3D<real>(nx, ny, 1, 0, "terrain_metric");
    metrics.terrain_slope_x =
        state::Field3D<real>(nx, ny, 1, 0, "terrain_slope_x_metric");
    metrics.terrain_slope_y =
        state::Field3D<real>(nx, ny, 1, 0, "terrain_slope_y_metric");
    metrics.z_interfaces_metric =
        state::Field3D<real>(nx, ny, nz + 1, 0, "z_interfaces_metric");
    metrics.z_centers_metric =
        state::Field3D<real>(nx, ny, nz, 0, "z_centers_metric");
    metrics.inv_dz_cell_metric =
        state::Field3D<real>(nx, ny, nz, 0, "inv_dz_cell_metric");
    metrics.inv_dz_face_metric =
        state::Field3D<real>(nx, ny, nz + 1, 0, "inv_dz_face_metric");

    const auto terrain_weight = [terrain_taper_eta](double eta) {
      if (eta <= 0.0) {
        return 1.0;
      }
      if (eta >= 1.0) {
        return 0.0;
      }
      if (terrain_taper_eta <= 0.0) {
        const double s = std::clamp(eta, 0.0, 1.0);
        return 1.0 - s * s * (3.0 - 2.0 * s);
      }
      if (eta <= terrain_taper_eta) {
        return 1.0;
      }
      const double s = std::clamp((eta - terrain_taper_eta) /
                                      std::max(1.0e-9, 1.0 - terrain_taper_eta),
                                  0.0, 1.0);
      return 1.0 - s * s * (3.0 - 2.0 * s);
    };

    for (int k = 0; k <= nz; ++k) {
      metrics.eta_interfaces[static_cast<std::size_t>(k)] =
          static_cast<double>(k) / static_cast<double>(nz);
      metrics.terrain_weights_interfaces[static_cast<std::size_t>(k)] =
          terrain_weight(metrics.eta_interfaces[static_cast<std::size_t>(k)]);
      metrics.z_interfaces[static_cast<std::size_t>(k)] =
          metrics.eta_interfaces[static_cast<std::size_t>(k)] * z_top;
      if (k < nz) {
        metrics.eta_centers[static_cast<std::size_t>(k)] =
            (static_cast<double>(k) + 0.5) / static_cast<double>(nz);
        metrics.terrain_weights_centers[static_cast<std::size_t>(k)] =
            terrain_weight(metrics.eta_centers[static_cast<std::size_t>(k)]);
        metrics.z_centers[static_cast<std::size_t>(k)] =
            metrics.eta_centers[static_cast<std::size_t>(k)] * z_top;
      }
    }

    auto terrain_at = [&](int ii, int jj) -> real {
      return terrain_dyn[static_cast<std::size_t>(jj) * nx +
                         static_cast<std::size_t>(ii)];
    };

    auto wrapped_i = [nx, periodic_x](int ii) {
      if (!periodic_x) {
        return std::clamp(ii, 0, nx - 1);
      }
      const int mod = ii % nx;
      return mod < 0 ? mod + nx : mod;
    };

    auto wrapped_j = [ny, periodic_y](int jj) {
      if (!periodic_y) {
        return std::clamp(jj, 0, ny - 1);
      }
      const int mod = jj % ny;
      return mod < 0 ? mod + ny : mod;
    };

    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        const real terrain_here = terrain_at(i, j);
        metrics.terrain(i, j, 0) = terrain_here;
        const int iw = wrapped_i(i - 1);
        const int ie = wrapped_i(i + 1);
        const int js = wrapped_j(j - 1);
        const int jn = wrapped_j(j + 1);

        const double ddx =
            (static_cast<double>(terrain_at(ie, j)) -
             static_cast<double>(terrain_at(iw, j))) /
            (periodic_x || (i > 0 && i + 1 < nx) ? (2.0 * dx) : dx);
        const double ddy =
            (static_cast<double>(terrain_at(i, jn)) -
             static_cast<double>(terrain_at(i, js))) /
            (periodic_y || (j > 0 && j + 1 < ny) ? (2.0 * dy) : dy);
        metrics.terrain_slope_x(i, j, 0) = static_cast<real>(ddx);
        metrics.terrain_slope_y(i, j, 0) = static_cast<real>(ddy);

        for (int k = 0; k <= nz; ++k) {
          metrics.z_interfaces_metric(i, j, k) = static_cast<real>(
              metrics.z_interfaces[static_cast<std::size_t>(k)] +
              metrics.terrain_weights_interfaces[static_cast<std::size_t>(k)] *
                  static_cast<double>(terrain_here));
          if (k == 0 || k == nz) {
            metrics.inv_dz_face_metric(i, j, k) = 0.0f;
          } else {
            const double dz_face =
                static_cast<double>(metrics.z_centers[static_cast<std::size_t>(k)]) +
                metrics.terrain_weights_centers[static_cast<std::size_t>(k)] *
                    static_cast<double>(terrain_here) -
                (static_cast<double>(metrics.z_centers[static_cast<std::size_t>(k - 1)]) +
                 metrics.terrain_weights_centers[static_cast<std::size_t>(k - 1)] *
                     static_cast<double>(terrain_here));
            metrics.inv_dz_face_metric(i, j, k) =
                static_cast<real>(1.0 / std::max(dz_face, 1.0e-6));
          }
        }

        for (int k = 0; k < nz; ++k) {
          metrics.z_centers_metric(i, j, k) = static_cast<real>(
              metrics.z_centers[static_cast<std::size_t>(k)] +
              metrics.terrain_weights_centers[static_cast<std::size_t>(k)] *
                  static_cast<double>(terrain_here));
          const double dz_cell =
              static_cast<double>(metrics.z_interfaces_metric(i, j, k + 1)) -
              static_cast<double>(metrics.z_interfaces_metric(i, j, k));
          metrics.inv_dz_cell_metric(i, j, k) =
              static_cast<real>(1.0 / std::max(dz_cell, 1.0e-6));
        }
      }
    }

    return metrics;
  }
};

}  // namespace gwm::domain
