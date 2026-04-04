#include <cstdlib>
#include <exception>
#include <string>

#include "gwm/core/moist_thermo.hpp"

#include "test_assert.hpp"

int main() {
  try {
    using namespace gwm;

    const real temperature_k = 300.0f;
    const real pressure_pa = 90000.0f;
    const real saturation_vapor_pressure =
        core::saturation_vapor_pressure_pa(temperature_k);
    const real expected_qsat =
        core::kEpsilon * saturation_vapor_pressure /
        (pressure_pa - (1.0f - core::kEpsilon) * saturation_vapor_pressure);
    const real qsat =
        core::saturation_specific_humidity(temperature_k, pressure_pa);

    TEST_NEAR(qsat, expected_qsat, 1.0e-6f);
    TEST_NEAR(
        core::relative_humidity_from_specific_humidity(qsat, temperature_k,
                                                       pressure_pa),
        100.0f, 1.0e-4f);
    TEST_NEAR(core::dewpoint_from_specific_humidity(qsat, pressure_pa),
              temperature_k, 2.0e-1f);

    return 0;
  } catch (const std::exception& ex) {
    test_fail(std::string("exception: ") + ex.what());
  }
}
