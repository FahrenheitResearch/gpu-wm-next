#include "gwm/state/tracer_registry.hpp"

#include "test_assert.hpp"

int main() {
  gwm::state::TracerRegistry registry;
  TEST_CHECK(registry.size() == 0);
  const int qv = registry.add({"qv", "kg/kg", true});
  const int qc = registry.add({"qc", "kg/kg", true});
  TEST_CHECK(qv == 0);
  TEST_CHECK(qc == 1);
  TEST_CHECK(registry.size() == 2);
  TEST_CHECK(registry.find("qv").value() == 0);
  TEST_CHECK(registry.find("qc").value() == 1);
  TEST_CHECK(!registry.find("qr").has_value());
  TEST_CHECK(registry.at(1).name == "qc");
  return 0;
}
