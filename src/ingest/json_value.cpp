#include "gwm/ingest/json_value.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

namespace gwm::ingest::json {

namespace {

Value from_nlohmann(const nlohmann::json& value) {
  if (value.is_null()) {
    return Value(nullptr);
  }
  if (value.is_boolean()) {
    return Value(value.get<bool>());
  }
  if (value.is_number()) {
    return Value(value.get<double>());
  }
  if (value.is_string()) {
    return Value(value.get<std::string>());
  }
  if (value.is_array()) {
    Value::Array array;
    array.reserve(value.size());
    for (const auto& element : value) {
      array.push_back(from_nlohmann(element));
    }
    return Value(std::move(array));
  }
  if (value.is_object()) {
    Value::Object object;
    object.reserve(value.size());
    for (const auto& [key, element] : value.items()) {
      object.emplace(key, from_nlohmann(element));
    }
    return Value(std::move(object));
  }
  throw std::runtime_error("Unsupported nlohmann::json value type");
}

}  // namespace

bool Value::as_bool() const {
  gwm::require(is_bool(), "JSON value is not bool");
  return std::get<bool>(storage_);
}

double Value::as_number() const {
  gwm::require(is_number(), "JSON value is not number");
  return std::get<double>(storage_);
}

int Value::as_int() const {
  const double value = as_number();
  gwm::require(std::fabs(value - std::round(value)) < 1.0e-6,
               "JSON numeric value is not integral");
  return static_cast<int>(std::llround(value));
}

std::uint64_t Value::as_uint64() const {
  const double value = as_number();
  gwm::require(value >= 0.0, "JSON numeric value is negative");
  gwm::require(std::fabs(value - std::round(value)) < 1.0e-6,
               "JSON numeric value is not integral");
  return static_cast<std::uint64_t>(std::llround(value));
}

const std::string& Value::as_string() const {
  gwm::require(is_string(), "JSON value is not string");
  return std::get<std::string>(storage_);
}

const Value::Array& Value::as_array() const {
  gwm::require(is_array(), "JSON value is not array");
  return std::get<Array>(storage_);
}

const Value::Object& Value::as_object() const {
  gwm::require(is_object(), "JSON value is not object");
  return std::get<Object>(storage_);
}

const Value& Value::at(const std::string& key) const {
  const auto& object = as_object();
  const auto it = object.find(key);
  gwm::require(it != object.end(), "Missing JSON object key: " + key);
  return it->second;
}

bool Value::has(const std::string& key) const {
  if (!is_object()) {
    return false;
  }
  return as_object().find(key) != as_object().end();
}

std::size_t Value::size() const {
  if (is_array()) {
    return as_array().size();
  }
  if (is_object()) {
    return as_object().size();
  }
  if (is_string()) {
    return as_string().size();
  }
  return 0;
}

Value parse_text(const std::string& text) {
  return from_nlohmann(nlohmann::json::parse(text));
}

Value parse_file(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("Failed to open JSON file: " + path);
  }
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return parse_text(buffer.str());
}

std::vector<real> to_real_vector(const Value& value) {
  const auto& array = value.as_array();
  std::vector<real> out;
  out.reserve(array.size());
  for (const auto& element : array) {
    out.push_back(static_cast<real>(element.as_number()));
  }
  return out;
}

std::unordered_map<std::string, std::string> to_string_map(const Value& value) {
  const auto& object = value.as_object();
  std::unordered_map<std::string, std::string> out;
  out.reserve(object.size());
  for (const auto& [key, entry] : object) {
    if (entry.is_string()) {
      out.emplace(key, entry.as_string());
    } else if (entry.is_number()) {
      out.emplace(key, std::to_string(entry.as_number()));
    } else if (entry.is_bool()) {
      out.emplace(key, entry.as_bool() ? "true" : "false");
    } else if (entry.is_null()) {
      out.emplace(key, "null");
    }
  }
  return out;
}

}  // namespace gwm::ingest::json
