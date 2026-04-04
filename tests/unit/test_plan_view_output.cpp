#include <algorithm>
#include <string>

#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/io/plan_view_output.hpp"

#include "test_assert.hpp"

namespace {

const gwm::io::PlanViewField& require_field(const gwm::io::PlanViewBundle& bundle,
                                            const std::string& name) {
  const auto it = std::find_if(
      bundle.fields.begin(), bundle.fields.end(),
      [&](const auto& field) { return field.name == name; });
  TEST_CHECK(it != bundle.fields.end());
  return *it;
}

}  // namespace

int main() {
  try {
    using namespace gwm;

    const auto layout =
        comm::VirtualRankLayout::build(6, 4, 3, 1, 2, 1, true, true);
    auto states = dycore::make_constant_dry_state(layout, 1.2f, 300.0f, 10.0f,
                                                  -4.0f, 0.5f,
                                                  "plan_view_test");
    const auto metrics = domain::GridMetrics::make_hybrid_height(
        6, 4, 3, 1000.0, 1000.0, 9000.0, std::vector<real>(24, 0.0f));

    const auto bundle = io::extract_dry_plan_view(states, layout, metrics,
                                                  "constant", 4, 2.0f, 1);
    TEST_CHECK(bundle.schema_version == "gwm-next-plan-view/v1");
    TEST_CHECK(bundle.case_kind == "constant");
    TEST_CHECK(bundle.steps == 4);
    TEST_NEAR(bundle.dt, 2.0f, 1.0e-6f);
    TEST_CHECK(bundle.slice_k == 1);
    TEST_CHECK(bundle.nx == 6);
    TEST_CHECK(bundle.ny == 4);
    TEST_CHECK(bundle.fields.size() >= 8);

    const auto& terrain = require_field(bundle, "terrain_height");
    const auto& rho = require_field(bundle, "rho_d");
    const auto& theta = require_field(bundle, "theta_m");
    const auto& u = require_field(bundle, "u_velocity");
    const auto& v = require_field(bundle, "v_velocity");
    const auto& wind = require_field(bundle, "wind_speed");
    const auto& w = require_field(bundle, "w_velocity");

    TEST_CHECK(static_cast<int>(terrain.values.size()) == bundle.nx * bundle.ny);
    TEST_CHECK(static_cast<int>(rho.values.size()) == bundle.nx * bundle.ny);
    TEST_NEAR(rho.values.front(), 1.2, 1.0e-6);
    TEST_NEAR(theta.values.front(), 300.0, 1.0e-6);
    TEST_NEAR(u.values.front(), 10.0, 1.0e-6);
    TEST_NEAR(v.values.front(), -4.0, 1.0e-6);
    TEST_NEAR(w.values.front(), 0.5, 1.0e-6);
    TEST_NEAR(wind.values.front(), std::sqrt(116.0), 1.0e-6);
    TEST_NEAR(terrain.values.front(), 0.0, 1.0e-6);

    const auto json = io::plan_view_bundle_to_json(bundle);
    TEST_CHECK(json.find("\"schema_version\": \"gwm-next-plan-view/v1\"") !=
               std::string::npos);
    TEST_CHECK(json.find("\"name\": \"theta_m\"") != std::string::npos);
    TEST_CHECK(json.find("\"storage\": \"row_major_yx\"") != std::string::npos);
    return 0;
  } catch (const std::exception& ex) {
    test_fail(std::string("exception: ") + ex.what());
  }
}
