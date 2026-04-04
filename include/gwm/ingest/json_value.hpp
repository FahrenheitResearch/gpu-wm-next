#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include "gwm/core/types.hpp"

namespace gwm::ingest::json {

class Value {
 public:
  using Object = std::unordered_map<std::string, Value>;
  using Array = std::vector<Value>;

  Value() = default;
  explicit Value(std::nullptr_t) : storage_(nullptr) {}
  explicit Value(bool value) : storage_(value) {}
  explicit Value(double value) : storage_(value) {}
  explicit Value(std::string value) : storage_(std::move(value)) {}
  explicit Value(Array value) : storage_(std::move(value)) {}
  explicit Value(Object value) : storage_(std::move(value)) {}

  [[nodiscard]] bool is_null() const {
    return std::holds_alternative<std::nullptr_t>(storage_);
  }
  [[nodiscard]] bool is_bool() const {
    return std::holds_alternative<bool>(storage_);
  }
  [[nodiscard]] bool is_number() const {
    return std::holds_alternative<double>(storage_);
  }
  [[nodiscard]] bool is_string() const {
    return std::holds_alternative<std::string>(storage_);
  }
  [[nodiscard]] bool is_array() const {
    return std::holds_alternative<Array>(storage_);
  }
  [[nodiscard]] bool is_object() const {
    return std::holds_alternative<Object>(storage_);
  }

  [[nodiscard]] bool as_bool() const;
  [[nodiscard]] double as_number() const;
  [[nodiscard]] int as_int() const;
  [[nodiscard]] std::uint64_t as_uint64() const;
  [[nodiscard]] const std::string& as_string() const;
  [[nodiscard]] const Array& as_array() const;
  [[nodiscard]] const Object& as_object() const;
  [[nodiscard]] const Value& at(const std::string& key) const;
  [[nodiscard]] bool has(const std::string& key) const;
  [[nodiscard]] std::size_t size() const;

 private:
  std::variant<std::nullptr_t, bool, double, std::string, Array, Object>
      storage_ = nullptr;
};

[[nodiscard]] Value parse_text(const std::string& text);
[[nodiscard]] Value parse_file(const std::string& path);
[[nodiscard]] std::vector<real> to_real_vector(const Value& value);
[[nodiscard]] std::unordered_map<std::string, std::string> to_string_map(
    const Value& value);

}  // namespace gwm::ingest::json
