#include <exception>

#include "gwm/state/tracer_registry.hpp"
#include "gwm/state/tracer_state.hpp"

#include "test_assert.hpp"

int main() {
  try {
    gwm::state::TracerRegistry registry;
    TEST_CHECK(registry.size() == 0);

    const int qv = registry.add({"specific_humidity", "kg/kg", true});
    const int qc = registry.add({"cloud_water_mixing_ratio", "kg/kg", true});
    const int signed_tracer = registry.add({"signed_tracer", "1", false});

    TEST_CHECK(qv == 0);
    TEST_CHECK(qc == 1);
    TEST_CHECK(signed_tracer == 2);
    TEST_CHECK(registry.size() == 3);

    TEST_CHECK(registry.find("specific_humidity").value() == 0);
    TEST_CHECK(registry.find("cloud_water_mixing_ratio").value() == 1);
    TEST_CHECK(registry.find("signed_tracer").value() == 2);
    TEST_CHECK(!registry.find("water_vapor_mixing_ratio").has_value());

    TEST_CHECK(registry.at(qv).name == "specific_humidity");
    TEST_CHECK(registry.at(qv).units == "kg/kg");
    TEST_CHECK(registry.at(qv).positive);
    TEST_CHECK(registry.at(qc).positive);
    TEST_CHECK(!registry.at(signed_tracer).positive);

    bool duplicate_rejected = false;
    try {
      static_cast<void>(
          registry.add({"specific_humidity", "kg/kg", true}));
    } catch (const std::exception&) {
      duplicate_rejected = true;
    }
    TEST_CHECK(duplicate_rejected);

    auto runtime_registry = gwm::state::make_specific_humidity_registry();
    gwm::state::TracerState tracer_state(std::move(runtime_registry), 3, 2, 2, 1,
                                         "runtime_tracer");
    TEST_CHECK(tracer_state.size() == 1);
    tracer_state.mass(gwm::state::kSpecificHumidityTracerName).fill(0.25f);
    TEST_NEAR(tracer_state.total_mass(0), 3.0, 1.0e-6);

    auto cloned = tracer_state.clone_empty_like("_clone");
    cloned.copy_all_from(tracer_state);
    TEST_NEAR(cloned.mass(0)(1, 1, 1), 0.25f, 1.0e-6f);

    auto warm_rain_registry = gwm::state::make_warm_rain_registry();
    TEST_CHECK(warm_rain_registry.size() == 3);
    TEST_CHECK(warm_rain_registry.find(gwm::state::kSpecificHumidityTracerName)
                   .value() == 0);
    TEST_CHECK(warm_rain_registry.find(gwm::state::kCloudWaterTracerName)
                   .value() == 1);
    TEST_CHECK(warm_rain_registry.find(gwm::state::kRainWaterTracerName)
                   .value() == 2);
    return 0;
  } catch (const std::exception& ex) {
    test_fail(std::string("exception: ") + ex.what());
  }
}
