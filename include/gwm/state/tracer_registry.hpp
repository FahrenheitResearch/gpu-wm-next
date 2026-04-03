#pragma once

#include <optional>
#include <string>
#include <vector>

#include "gwm/core/types.hpp"

namespace gwm::state {

struct TracerSpec {
  std::string name;
  std::string units;
  bool positive = true;
};

class TracerRegistry {
 public:
  int add(TracerSpec spec) {
    gwm::require(find(spec.name) == std::nullopt,
                 "Tracer already exists in registry");
    specs_.push_back(std::move(spec));
    return static_cast<int>(specs_.size() - 1);
  }

  [[nodiscard]] std::optional<int> find(const std::string& name) const {
    for (std::size_t i = 0; i < specs_.size(); ++i) {
      if (specs_[i].name == name) {
        return static_cast<int>(i);
      }
    }
    return std::nullopt;
  }

  [[nodiscard]] const TracerSpec& at(int index) const {
    return specs_.at(static_cast<std::size_t>(index));
  }

  [[nodiscard]] std::size_t size() const { return specs_.size(); }

 private:
  std::vector<TracerSpec> specs_;
};

}  // namespace gwm::state
