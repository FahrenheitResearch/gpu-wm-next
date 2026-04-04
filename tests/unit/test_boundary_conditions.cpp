#include "gwm/comm/virtual_rank_layout.hpp"
#include "gwm/dycore/boundary_conditions.hpp"

#include "test_assert.hpp"

int main() {
  using namespace gwm;

  const auto layout =
      comm::VirtualRankLayout::build(5, 4, 3, 1, 1, 1, false, false);
  TEST_CHECK(layout.size() == 1);
  TEST_CHECK(dycore::touches_west_boundary(layout.front()));
  TEST_CHECK(dycore::touches_east_boundary(layout.front()));
  TEST_CHECK(dycore::touches_south_boundary(layout.front()));
  TEST_CHECK(dycore::touches_north_boundary(layout.front()));

  auto states = dycore::make_constant_dry_state(layout, 1.0f, 300.0f, 0.0f,
                                                0.0f, 0.0f, "state");
  auto boundary_states = dycore::make_constant_dry_state(
      layout, 1.5f, 305.0f, 4.0f, -3.0f, 0.5f, "boundary");

  std::vector<state::TracerState> tracer_states;
  std::vector<state::TracerState> tracer_boundary_states;
  tracer_states.emplace_back(state::make_specific_humidity_registry(),
                             layout.front().nx_local(),
                             layout.front().ny_local(), layout.front().nz,
                             layout.front().halo, "tracer_state");
  tracer_boundary_states.emplace_back(state::make_specific_humidity_registry(),
                                      layout.front().nx_local(),
                                      layout.front().ny_local(),
                                      layout.front().nz, layout.front().halo,
                                      "tracer_boundary");
  tracer_states.front().mass(gwm::state::kSpecificHumidityTracerName).fill(
      0.01f);
  tracer_boundary_states.front().mass(gwm::state::kSpecificHumidityTracerName)
      .fill(0.015f);

  dycore::apply_reference_boundaries(states, boundary_states, layout);
  dycore::apply_reference_boundaries(tracer_states, tracer_boundary_states,
                                     layout);

  auto& state = states.front();
  auto& tracer_state = tracer_states.front();

  TEST_NEAR(state.rho_d(0, 0, 0), 1.5f, 1.0e-6f);
  TEST_NEAR(state.rho_theta_m(0, 0, 0), 1.5f * 305.0f, 1.0e-6f);
  TEST_NEAR(state.mom_u.storage()(0, 0, 0), 1.5f * 4.0f, 1.0e-6f);
  TEST_NEAR(state.mom_v.storage()(0, 0, 0), 1.5f * -3.0f, 1.0e-6f);
  TEST_NEAR(state.mom_w.storage()(0, 0, 0), 1.5f * 0.5f, 1.0e-6f);

  TEST_NEAR(state.rho_d(state.rho_d.nx() - 1, state.rho_d.ny() - 1,
                        state.rho_d.nz() - 1),
            1.5f, 1.0e-6f);
  TEST_NEAR(state.mom_u.storage()(state.mom_u.storage().nx() - 1,
                                  state.mom_u.storage().ny() - 1,
                                  state.mom_u.storage().nz() - 1),
            1.5f * 4.0f, 1.0e-6f);

  TEST_NEAR(state.rho_d(2, 2, 1), 1.0f, 1.0e-6f);
  TEST_NEAR(state.rho_theta_m(2, 2, 1), 300.0f, 1.0e-6f);
  TEST_NEAR(state.mom_u.storage()(2, 2, 1), 0.0f, 1.0e-6f);
  TEST_NEAR(state.mom_v.storage()(2, 2, 1), 0.0f, 1.0e-6f);
  TEST_NEAR(state.mom_w.storage()(2, 2, 1), 0.0f, 1.0e-6f);
  TEST_NEAR(
      tracer_state.mass(gwm::state::kSpecificHumidityTracerName)(0, 0, 0),
      0.015f,
      1.0e-6f);
  TEST_NEAR(
      tracer_state.mass(gwm::state::kSpecificHumidityTracerName)(2, 2, 1),
      0.01f,
      1.0e-6f);
  return 0;
}
