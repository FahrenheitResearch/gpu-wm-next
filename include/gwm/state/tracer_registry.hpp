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

inline constexpr const char* kSpecificHumidityTracerName = "specific_humidity";
inline constexpr const char* kCloudWaterTracerName =
    "cloud_water_mixing_ratio";
inline constexpr const char* kRainWaterTracerName =
    "rain_water_mixing_ratio";

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

[[nodiscard]] inline TracerRegistry make_specific_humidity_registry() {
  TracerRegistry registry;
  registry.add({kSpecificHumidityTracerName, "kg/kg", true});
  return registry;
}

[[nodiscard]] inline TracerRegistry make_warm_rain_registry() {
  TracerRegistry registry;
  registry.add({kSpecificHumidityTracerName, "kg/kg", true});
  registry.add({kCloudWaterTracerName, "kg/kg", true});
  registry.add({kRainWaterTracerName, "kg/kg", true});
  return registry;
}

}  // namespace gwm::state
