#pragma once

#include <cuda_runtime.h>

#include <sstream>
#include <stdexcept>

namespace gwm::core {

inline void check_cuda(cudaError_t status, const char* expr, const char* file,
                       int line) {
  if (status == cudaSuccess) {
    return;
  }
  std::ostringstream oss;
  oss << "CUDA error at " << file << ":" << line << " for " << expr << ": "
      << cudaGetErrorString(status);
  throw std::runtime_error(oss.str());
}

}  // namespace gwm::core

#define GWM_CUDA_CHECK(expr) \
  ::gwm::core::check_cuda((expr), #expr, __FILE__, __LINE__)
