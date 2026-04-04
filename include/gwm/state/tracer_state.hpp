#pragma once

#include <optional>
#include <string>
#include <vector>

#include "gwm/state/field3d.hpp"
#include "gwm/state/tracer_registry.hpp"

namespace gwm::state {

class TracerState {
 public:
  TracerState() = default;

  TracerState(TracerRegistry registry, int nx, int ny, int nz, int halo,
              std::string label_prefix = {})
      : registry_(std::move(registry)),
        label_prefix_(label_prefix.empty() ? "tracer_state"
                                           : std::move(label_prefix)) {
    masses_.reserve(registry_.size());
    for (std::size_t idx = 0; idx < registry_.size(); ++idx) {
      masses_.emplace_back(nx, ny, nz, halo,
                           label_prefix_ + "_" + registry_.at(
                               static_cast<int>(idx)).name);
    }
  }

  [[nodiscard]] int nx() const {
    return masses_.empty() ? 0 : masses_.front().nx();
  }
  [[nodiscard]] int ny() const {
    return masses_.empty() ? 0 : masses_.front().ny();
  }
  [[nodiscard]] int nz() const {
    return masses_.empty() ? 0 : masses_.front().nz();
  }
  [[nodiscard]] int halo() const {
    return masses_.empty() ? 0 : masses_.front().halo();
  }

  [[nodiscard]] const TracerRegistry& registry() const { return registry_; }
  [[nodiscard]] std::size_t size() const { return masses_.size(); }

  [[nodiscard]] std::optional<int> find(const std::string& name) const {
    return registry_.find(name);
  }

  [[nodiscard]] Field3D<real>& mass(int index) {
    return masses_.at(static_cast<std::size_t>(index));
  }

  [[nodiscard]] const Field3D<real>& mass(int index) const {
    return masses_.at(static_cast<std::size_t>(index));
  }

  [[nodiscard]] Field3D<real>& mass(const std::string& name) {
    const auto index = registry_.find(name);
    gwm::require(index.has_value(), "Tracer not found in TracerState: " + name);
    return mass(*index);
  }

  [[nodiscard]] const Field3D<real>& mass(const std::string& name) const {
    const auto index = registry_.find(name);
    gwm::require(index.has_value(), "Tracer not found in TracerState: " + name);
    return mass(*index);
  }

  [[nodiscard]] TracerState clone_empty_like(
      const std::string& suffix = "_tmp") const {
    return TracerState(registry_, nx(), ny(), nz(), halo(),
                       label_prefix_ + suffix);
  }

  void copy_all_from(const TracerState& other) {
    require_compatible(other, "copy_all_from");
    for (std::size_t idx = 0; idx < masses_.size(); ++idx) {
      masses_[idx].copy_all_from(other.masses_[idx]);
    }
  }

  void fill_zero() {
    for (auto& field : masses_) {
      field.fill(0.0f);
    }
  }

  [[nodiscard]] double total_mass(int index) const {
    return mass(index).owned_sum();
  }

 private:
  TracerRegistry registry_{};
  std::vector<Field3D<real>> masses_{};
  std::string label_prefix_ = "tracer_state";

  void require_compatible(const TracerState& other,
                          const char* context) const {
    gwm::require(registry_.size() == other.registry_.size(),
                 std::string("Tracer count mismatch in TracerState::") +
                     context);
    gwm::require(masses_.size() == other.masses_.size(),
                 std::string("Tracer storage mismatch in TracerState::") +
                     context);
    for (std::size_t idx = 0; idx < registry_.size(); ++idx) {
      const auto spec_index = static_cast<int>(idx);
      gwm::require(registry_.at(spec_index).name ==
                       other.registry_.at(spec_index).name,
                   std::string("Tracer name mismatch in TracerState::") +
                       context);
      gwm::require(masses_[idx].nx() == other.masses_[idx].nx() &&
                       masses_[idx].ny() == other.masses_[idx].ny() &&
                       masses_[idx].nz() == other.masses_[idx].nz() &&
                       masses_[idx].halo() == other.masses_[idx].halo(),
                   std::string("Tracer field shape mismatch in TracerState::") +
                       context);
    }
  }
};

}  // namespace gwm::state
