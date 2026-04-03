#pragma once

#include <string>

namespace gwm::comm {

struct MpiRuntimeInfo {
  bool compiled_with_mpi = false;
  bool launcher_available = false;
  std::string backend_name;
};

[[nodiscard]] MpiRuntimeInfo query_mpi_runtime();

}  // namespace gwm::comm
