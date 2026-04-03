#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace gwm {

using real = float;
using index_t = int;

struct Int2 {
  int x = 0;
  int y = 0;
};

struct Int3 {
  int x = 0;
  int y = 0;
  int z = 0;
};

inline void require(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

}  // namespace gwm
