#pragma once

#include <memory>
#include <string>
#include <vector>

#include "gwm/ingest/canonical_ir.hpp"

namespace gwm::ingest {

struct ExternalToolDescriptor {
  std::string name;
  std::string local_repo_hint;
  std::string role;
  bool allowed_in_runtime_core = false;
};

[[nodiscard]] std::string to_string(SourceKind kind);
[[nodiscard]] SourceKind source_kind_from_string(const std::string& value);
[[nodiscard]] bool is_first_class_source(SourceKind kind);
[[nodiscard]] std::unique_ptr<SourceAdapter> make_adapter(SourceKind kind);
[[nodiscard]] std::vector<std::unique_ptr<SourceAdapter>> make_default_adapters();
[[nodiscard]] std::vector<ExternalToolDescriptor> recommended_external_tools();

}  // namespace gwm::ingest
