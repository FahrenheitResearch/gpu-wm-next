#include "gwm/ingest/source_catalog.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm::ingest;

  TEST_CHECK(to_string(SourceKind::HRRR) == "HRRR");
  TEST_CHECK(to_string(SourceKind::RRFS) == "RRFS");

  auto hrrr = make_adapter(SourceKind::HRRR);
  TEST_CHECK(hrrr != nullptr);
  TEST_CHECK(hrrr->name() == "HRRR");
  TEST_CHECK(hrrr->boundary_interval_seconds() == 3600);
  TEST_CHECK(!hrrr->atmosphere_field_bindings().empty());
  TEST_CHECK(!hrrr->surface_field_bindings().empty());
  TEST_CHECK(!hrrr->static_field_bindings().empty());

  auto rrfs = make_adapter(SourceKind::RRFS);
  TEST_CHECK(rrfs != nullptr);
  TEST_CHECK(rrfs->name() == "RRFS");
  TEST_CHECK(rrfs->boundary_interval_seconds() == 3600);
  TEST_CHECK(!rrfs->atmosphere_field_bindings().empty());

  auto gfs = make_adapter(SourceKind::GFS);
  TEST_CHECK(gfs == nullptr);
  TEST_CHECK(is_first_class_source(SourceKind::HRRR));
  TEST_CHECK(is_first_class_source(SourceKind::RRFS));
  TEST_CHECK(!is_first_class_source(SourceKind::GFS));

  const auto adapters = make_default_adapters();
  TEST_CHECK(adapters.size() == 2);

  const auto tools = recommended_external_tools();
  TEST_CHECK(!tools.empty());
  for (const auto& tool : tools) {
    TEST_CHECK(!tool.name.empty());
    TEST_CHECK(!tool.role.empty());
    TEST_CHECK(tool.allowed_in_runtime_core == false);
  }

  return 0;
}
