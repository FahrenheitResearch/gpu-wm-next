#include "gwm/surface/surface_layer_exchange.hpp"

#include <algorithm>

namespace gwm::surface {

namespace {

real tile_surface_humidity(const SurfaceState& state,
                           const SurfaceLayerExchangeForcing& forcing, int tile,
                           int i, int j) {
  if (!forcing.tile_q_surface.empty()) {
    gwm::require(
        static_cast<int>(forcing.tile_q_surface.size()) == state.ntile(),
        "Surface-layer exchange tile_q_surface must be empty or sized to ntile");
    return std::max(forcing.tile_q_surface.at(static_cast<std::size_t>(tile)),
                    0.0f);
  }

  return saturation_specific_humidity(state.skin_temperature(tile, i, j),
                                      forcing.psfc);
}

void accumulate_weighted(SurfaceLayerDiagnostics& cell_mean,
                         const SurfaceLayerDiagnostics& tile_diag, real weight) {
  cell_mean.t2 += weight * tile_diag.t2;
  cell_mean.q2 += weight * tile_diag.q2;
  cell_mean.rh2 += weight * tile_diag.rh2;
  cell_mean.u10 += weight * tile_diag.u10;
  cell_mean.v10 += weight * tile_diag.v10;
  cell_mean.wind_speed_ref += weight * tile_diag.wind_speed_ref;
  cell_mean.wind_speed_target += weight * tile_diag.wind_speed_target;
  cell_mean.cm += weight * tile_diag.cm;
  cell_mean.ch += weight * tile_diag.ch;
  cell_mean.cq += weight * tile_diag.cq;
  cell_mean.ustar += weight * tile_diag.ustar;
}

}  // namespace

SurfaceLayerExchangeResult evaluate_neutral_surface_exchange(
    const SurfaceState& state, const SurfaceStaticProperties& properties, int i,
    int j, const SurfaceLayerExchangeForcing& forcing) {
  gwm::require(state.nx() == properties.nx() && state.ny() == properties.ny() &&
                   state.ntile() == properties.ntile(),
               "Surface state and static properties must share nx/ny/ntile");
  gwm::require(state.ntile() > 0, "Surface-layer exchange requires ntile > 0");
  properties.validate_cell_tile_fractions(i, j);

  SurfaceLayerExchangeResult out{};
  out.tile_outputs.reserve(static_cast<std::size_t>(state.ntile()));

  for (int tile = 0; tile < state.ntile(); ++tile) {
    SurfaceLayerInputs tile_inputs{};
    tile_inputs.z_ref = forcing.z_ref;
    tile_inputs.z_target_temp = forcing.z_target_temp;
    tile_inputs.z_target_wind = forcing.z_target_wind;
    tile_inputs.z0m = properties.z0m(tile, i, j);
    tile_inputs.z0h = properties.z0h(tile, i, j);
    tile_inputs.u_ref = forcing.u_ref;
    tile_inputs.v_ref = forcing.v_ref;
    tile_inputs.theta_ref = forcing.theta_ref;
    tile_inputs.q_ref = forcing.q_ref;
    tile_inputs.tskin = state.skin_temperature(tile, i, j);
    tile_inputs.q_surface = tile_surface_humidity(state, forcing, tile, i, j);
    tile_inputs.psfc = forcing.psfc;

    const auto diag = evaluate_neutral_surface_layer(tile_inputs);
    out.tile_outputs.push_back(diag);
    accumulate_weighted(out.cell_mean, diag, properties.tile_fraction(tile, i, j));
  }

  return out;
}

}  // namespace gwm::surface
