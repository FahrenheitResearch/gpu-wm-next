#pragma once

#include "gwm/core/types.hpp"

namespace gwm::domain {

struct SubdomainDescriptor {
  int rank = 0;
  int coord_x = 0;
  int coord_y = 0;
  int ranks_x = 1;
  int ranks_y = 1;
  int global_nx = 0;
  int global_ny = 0;
  int nz = 0;
  int halo = 0;
  int i_begin = 0;
  int i_end = 0;
  int j_begin = 0;
  int j_end = 0;
  bool periodic_x = true;
  bool periodic_y = true;

  [[nodiscard]] int nx_local() const { return i_end - i_begin; }
  [[nodiscard]] int ny_local() const { return j_end - j_begin; }
};

}  // namespace gwm::domain
