#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

#include "gwm/core/cuda_utils.hpp"
#include "gwm/core/runtime_summary.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/dry_diagnostics.hpp"
#include "gwm/dycore/passive_tracer.hpp"
#include "gwm/ingest/prepared_case_init.hpp"
#include "gwm/ingest/runtime_case.hpp"
#include "gwm/ingest/source_catalog.hpp"
#include "gwm/io/plan_view_output.hpp"
#include "gwm/physics/warm_rain.hpp"
#include "gwm/surface/surface_runtime_init.hpp"

namespace {

struct DriverOptions {
  std::string analysis_state_path;
  std::string boundary_cache_path;
  gwm::ingest::PreparedCaseInitConfig init_config{};
  gwm::dycore::DryStepperConfig step_config{};
  gwm::physics::WarmRainConfig warm_rain_config = [] {
    gwm::physics::WarmRainConfig config{};
    config.rain_terminal_velocity = 12.0f;
    return config;
  }();
  int steps = 1;
  bool enable_boundaries = true;
  bool enable_fast_modes = true;
  bool enable_tracer_transport = true;
  bool enable_warm_rain = true;
  std::string summary_json_path;
  std::string plan_view_json_path;
  int plan_view_level = 0;
};

class NullFastModeIntegrator final : public gwm::dycore::FastModeIntegrator {
 public:
  void apply(std::vector<gwm::dycore::DryState>&,
             const std::vector<gwm::domain::SubdomainDescriptor>&,
             const gwm::domain::GridMetrics&,
             const gwm::dycore::DryStepperConfig&) override {}
};

bool driver_trace_enabled() {
  const char* value = std::getenv("GWM_TRACE_DRIVER");
  if (value == nullptr) {
    return false;
  }
  return std::string(value) != "0" && std::string(value) != "false";
}

std::optional<int> parse_optional_env_int(const char* name) {
  const char* value = std::getenv(name);
  if (value == nullptr || *value == '\0') {
    return std::nullopt;
  }
  return std::stoi(value);
}

std::string phase_summary_dir() {
  const char* value = std::getenv("GWM_PHASE_SUMMARY_DIR");
  return value == nullptr ? std::string{} : std::string(value);
}

bool phase_summary_step_selected(int step) {
  const auto exact_step = parse_optional_env_int("GWM_PHASE_SUMMARY_STEP");
  if (exact_step.has_value()) {
    return step == *exact_step;
  }
  const auto min_step = parse_optional_env_int("GWM_PHASE_SUMMARY_STEP_MIN");
  const auto max_step = parse_optional_env_int("GWM_PHASE_SUMMARY_STEP_MAX");
  if (min_step.has_value() && step < *min_step) {
    return false;
  }
  if (max_step.has_value() && step > *max_step) {
    return false;
  }
  return true;
}

void trace_driver(const std::string& message) {
  if (!driver_trace_enabled()) {
    return;
  }
  std::cerr << "[gwm_prepared_case_driver] " << message << std::endl;
}

void maybe_write_phase_summary(
    int step, gwm::real step_start_seconds, const std::string& phase,
    const std::vector<gwm::dycore::DryState>& states,
    const std::vector<gwm::state::TracerState>& tracers,
    const std::vector<gwm::domain::SubdomainDescriptor>& layout,
    const gwm::domain::GridMetrics& metrics,
    const std::vector<gwm::physics::WarmRainSurfaceAccumulation>&
        surface_precip_accum) {
  const auto dir = phase_summary_dir();
  if (dir.empty() || !phase_summary_step_selected(step)) {
    return;
  }

  std::filesystem::create_directories(dir);
  const auto summary = gwm::core::summarize_runtime_state(
      states, tracers, layout, metrics, &surface_precip_accum);
  std::ostringstream json;
  json << "{\n";
  json << "  \"step\": " << step << ",\n";
  json << "  \"step_start_seconds\": " << step_start_seconds << ",\n";
  json << "  \"phase\": \"" << phase << "\",\n";
  json << "  \"summary\": "
       << gwm::core::runtime_state_summary_to_json(summary, "    ") << "\n";
  json << "}\n";

  std::ostringstream filename;
  filename << "step_" << step << "_" << phase << ".json";
  const auto path = std::filesystem::path(dir) / filename.str();
  std::ofstream out(path, std::ios::binary);
  if (!out) {
    throw std::runtime_error("Failed to open phase summary path: " +
                             path.string());
  }
  out << json.str();
}

template <typename T>
T parse_value(const std::string& value);

template <>
int parse_value<int>(const std::string& value) {
  return std::stoi(value);
}

template <>
float parse_value<float>(const std::string& value) {
  return std::stof(value);
}

template <>
double parse_value<double>(const std::string& value) {
  return std::stod(value);
}

bool parse_bool(const std::string& value) {
  if (value == "true" || value == "1") {
    return true;
  }
  if (value == "false" || value == "0") {
    return false;
  }
  throw std::runtime_error("Invalid boolean value: " + value);
}

std::string require_value(int& index, int argc, char** argv) {
  if (index + 1 >= argc) {
    throw std::runtime_error(std::string("Missing value for ") + argv[index]);
  }
  ++index;
  return argv[index];
}

void parse_args(int argc, char** argv, DriverOptions& opts) {
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--analysis-state") {
      opts.analysis_state_path = require_value(i, argc, argv);
    } else if (arg == "--boundary-cache") {
      opts.boundary_cache_path = require_value(i, argc, argv);
    } else if (arg == "--steps") {
      opts.steps = parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--enable-boundaries") {
      opts.enable_boundaries = parse_bool(require_value(i, argc, argv));
    } else if (arg == "--enable-fast-modes") {
      opts.enable_fast_modes = parse_bool(require_value(i, argc, argv));
    } else if (arg == "--enable-tracer-transport") {
      opts.enable_tracer_transport = parse_bool(require_value(i, argc, argv));
    } else if (arg == "--enable-warm-rain") {
      opts.enable_warm_rain = parse_bool(require_value(i, argc, argv));
    } else if (arg == "--warm-rain-condensation-relaxation") {
      opts.warm_rain_config.condensation_relaxation =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--warm-rain-cloud-autoconversion-threshold") {
      opts.warm_rain_config.cloud_autoconversion_threshold =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--warm-rain-cloud-autoconversion-rate") {
      opts.warm_rain_config.cloud_autoconversion_rate =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--warm-rain-rain-evaporation-rate") {
      opts.warm_rain_config.rain_evaporation_rate =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--warm-rain-terminal-velocity") {
      opts.warm_rain_config.rain_terminal_velocity =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--warm-rain-enable-latent-heating") {
      opts.warm_rain_config.enable_latent_heating =
          parse_bool(require_value(i, argc, argv));
    } else if (arg == "--dt") {
      opts.step_config.dt = parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--fast-substeps") {
      opts.step_config.fast_substeps =
          parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--gravity") {
      opts.step_config.gravity = parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--fast-update-horizontal-momentum") {
      opts.step_config.fast_update_horizontal_momentum =
          parse_bool(require_value(i, argc, argv));
    } else if (arg == "--fast-update-density") {
      opts.step_config.fast_update_density =
          parse_bool(require_value(i, argc, argv));
    } else if (arg == "--halo") {
      opts.init_config.halo = parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--ranks-x") {
      opts.init_config.ranks_x = parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--ranks-y") {
      opts.init_config.ranks_y = parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--terrain-taper-eta") {
      opts.init_config.terrain_taper_eta =
          parse_value<double>(require_value(i, argc, argv));
    } else if (arg == "--periodic-x") {
      opts.init_config.periodic_x = parse_bool(require_value(i, argc, argv));
    } else if (arg == "--periodic-y") {
      opts.init_config.periodic_y = parse_bool(require_value(i, argc, argv));
    } else if (arg == "--summary-json") {
      opts.summary_json_path = require_value(i, argc, argv);
    } else if (arg == "--plan-view-json") {
      opts.plan_view_json_path = require_value(i, argc, argv);
    } else if (arg == "--plan-view-level") {
      opts.plan_view_level = parse_value<int>(require_value(i, argc, argv));
    } else {
      throw std::runtime_error("Unknown argument: " + arg);
    }
  }

  if (opts.analysis_state_path.empty()) {
    throw std::runtime_error("--analysis-state is required");
  }
  if (opts.boundary_cache_path.empty()) {
    throw std::runtime_error("--boundary-cache is required");
  }
  if (opts.steps < 0) {
    throw std::runtime_error("--steps must be nonnegative");
  }
}

std::string prepared_summary_json(
    const gwm::ingest::PreparedRuntimeCase& runtime_case, int steps,
    gwm::real dt, int fast_substeps, bool enable_boundaries,
    bool enable_fast_modes, bool enable_tracer_transport,
    bool enable_warm_rain, bool fast_update_horizontal_momentum,
    bool fast_update_density,
    const gwm::physics::WarmRainConfig& warm_rain_config,
    const gwm::ingest::PreparedCaseBalanceDiagnostics& startup_balance,
    const gwm::surface::SurfaceRuntimeInitResult& surface_runtime,
    const gwm::core::RuntimeStateSummary& initial,
    const gwm::core::RuntimeStateSummary& final) {
  const double elapsed_seconds =
      static_cast<double>(steps) * static_cast<double>(dt);
  const double elapsed_hours = elapsed_seconds / 3600.0;
  std::ostringstream oss;
  oss << std::setprecision(10);
  oss << "{\n";
  oss << "  \"case\": \"prepared_case\",\n";
  oss << "  \"source\": \"" << gwm::ingest::to_string(runtime_case.analysis.source)
      << "\",\n";
  oss << "  \"cycle_time_utc\": \"" << runtime_case.analysis.cycle_time_utc
      << "\",\n";
  oss << "  \"analysis_valid_time_utc\": \"" << runtime_case.analysis.valid_time_utc
      << "\",\n";
  oss << "  \"steps\": " << steps << ",\n";
  oss << "  \"dt\": " << dt << ",\n";
  oss << "  \"fast_substeps\": " << fast_substeps << ",\n";
  oss << "  \"enable_boundaries\": " << (enable_boundaries ? "true" : "false")
      << ",\n";
  oss << "  \"enable_fast_modes\": " << (enable_fast_modes ? "true" : "false")
      << ",\n";
  oss << "  \"enable_tracer_transport\": "
      << (enable_tracer_transport ? "true" : "false") << ",\n";
  oss << "  \"enable_warm_rain\": " << (enable_warm_rain ? "true" : "false")
      << ",\n";
  oss << "  \"fast_update_horizontal_momentum\": "
      << (fast_update_horizontal_momentum ? "true" : "false") << ",\n";
  oss << "  \"fast_update_density\": "
      << (fast_update_density ? "true" : "false") << ",\n";
  oss << "  \"warm_rain\": {\n";
  oss << "    \"condensation_relaxation\": "
      << warm_rain_config.condensation_relaxation << ",\n";
  oss << "    \"cloud_autoconversion_threshold\": "
      << warm_rain_config.cloud_autoconversion_threshold << ",\n";
  oss << "    \"cloud_autoconversion_rate\": "
      << warm_rain_config.cloud_autoconversion_rate << ",\n";
  oss << "    \"rain_evaporation_rate\": "
      << warm_rain_config.rain_evaporation_rate << ",\n";
  oss << "    \"rain_terminal_velocity\": "
      << warm_rain_config.rain_terminal_velocity << ",\n";
  oss << "    \"sedimentation_layer_depth_m\": "
      << warm_rain_config.sedimentation_layer_depth_m << ",\n";
  oss << "    \"enable_latent_heating\": "
      << (warm_rain_config.enable_latent_heating ? "true" : "false") << "\n";
  oss << "  },\n";
  oss << "  \"elapsed_seconds\": " << elapsed_seconds << ",\n";
  oss << "  \"elapsed_hours\": " << elapsed_hours << ",\n";
  oss << "  \"grid\": {\n";
  oss << "    \"nx\": " << runtime_case.analysis.grid.nx << ",\n";
  oss << "    \"ny\": " << runtime_case.analysis.grid.ny << ",\n";
  oss << "    \"nz\": " << runtime_case.analysis.grid.nz << ",\n";
  oss << "    \"dx\": " << runtime_case.analysis.grid.dx << ",\n";
  oss << "    \"dy\": " << runtime_case.analysis.grid.dy << ",\n";
  oss << "    \"z_top\": " << runtime_case.analysis.grid.z_top << "\n";
  oss << "  },\n";
  oss << "  \"surface_runtime\": {\n";
  oss << "    \"ntile\": " << surface_runtime.state.ntile() << ",\n";
  oss << "    \"nsoil\": " << surface_runtime.state.nsoil() << "\n";
  oss << "  },\n";
  oss << "  \"startup_balance\": {\n";
  oss << "    \"max_rel_eos\": " << startup_balance.max_rel_eos << ",\n";
  oss << "    \"max_rel_hydrostatic\": " << startup_balance.max_rel_hydrostatic
      << ",\n";
  oss << "    \"max_rel_fast_vertical\": "
      << startup_balance.max_rel_fast_vertical << ",\n";
  oss << "    \"max_abs_mass_divergence\": "
      << startup_balance.max_abs_mass_divergence << ",\n";
  oss << "    \"max_rel_mass_divergence\": "
      << startup_balance.max_rel_mass_divergence << ",\n";
  oss << "    \"max_abs_mom_w_bottom\": "
      << startup_balance.max_abs_mom_w_bottom << ",\n";
  oss << "    \"max_abs_mom_w_top\": " << startup_balance.max_abs_mom_w_top
      << ",\n";
  oss << "    \"max_abs_tracer_closure\": "
      << startup_balance.max_abs_tracer_closure;
  if (startup_balance.max_abs_z_src_minus_metric.has_value()) {
    oss << ",\n    \"max_abs_z_src_minus_metric\": "
        << *startup_balance.max_abs_z_src_minus_metric << "\n";
  } else {
    oss << "\n";
  }
  oss << "  },\n";
  oss << "  \"initial\": "
      << gwm::core::runtime_state_summary_to_json(initial, "    ") << ",\n";
  oss << "  \"final\": "
      << gwm::core::runtime_state_summary_to_json(final, "    ") << "\n";
  oss << "}";
  return oss.str();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    trace_driver("startup");
    DriverOptions opts{};
    parse_args(argc, argv, opts);
    trace_driver("args_parsed");

    const auto runtime_case = gwm::ingest::load_prepared_runtime_case(
        opts.analysis_state_path, opts.boundary_cache_path);
    trace_driver("runtime_case_loaded");
    const auto layout =
        gwm::ingest::make_prepared_case_layout(runtime_case.analysis,
                                               opts.init_config);
    trace_driver("layout_built");
    const auto metrics =
        gwm::ingest::make_prepared_case_metrics(runtime_case.analysis,
                                                opts.init_config);
    trace_driver("metrics_built");
    auto states = gwm::ingest::make_dry_states_from_analysis(
        runtime_case.analysis, metrics, layout, "prepared_case");
    trace_driver("dry_states_initialized");
    auto tracers = gwm::ingest::make_warm_rain_tracers_from_analysis(
        runtime_case.analysis, states, layout, "prepared_case_tracer");
    trace_driver("warm_rain_tracers_initialized");
    const auto startup_balance = gwm::ingest::diagnose_prepared_case_balance(
        runtime_case.analysis, metrics, states, layout, &tracers);
    trace_driver("startup_balance_done");

    const auto surface_runtime = gwm::surface::make_surface_runtime_from_canonical_fields(
        runtime_case.analysis.surface, runtime_case.analysis.static_surface,
        runtime_case.analysis.grid.nx, runtime_case.analysis.grid.ny);
    trace_driver("surface_runtime_initialized");
    (void)surface_runtime.properties;
    std::vector<gwm::physics::WarmRainSurfaceAccumulation> surface_precip_accum;
    surface_precip_accum.resize(states.size());
    for (std::size_t rank = 0; rank < states.size(); ++rank) {
      surface_precip_accum[rank].reset(states[rank].rho_d.nx(),
                                       states[rank].rho_d.ny(),
                                       metrics.dx * metrics.dy);
    }

    const auto initial =
        gwm::core::summarize_runtime_state(states, tracers, layout, metrics,
                                          &surface_precip_accum);
    trace_driver("initial_summary_done");

    gwm::ingest::PreparedCaseBoundaryUpdater prepared_boundary_updater(
        runtime_case.analysis, runtime_case.boundary_cache, opts.init_config);
    gwm::dycore::NullBoundaryUpdater null_boundary_updater;
    gwm::dycore::BoundaryUpdater& boundary_updater =
        opts.enable_boundaries
            ? static_cast<gwm::dycore::BoundaryUpdater&>(
                  prepared_boundary_updater)
            : static_cast<gwm::dycore::BoundaryUpdater&>(null_boundary_updater);
    gwm::ingest::PreparedCaseTracerBoundaryUpdater prepared_tracer_boundary_updater(
        runtime_case.analysis, runtime_case.boundary_cache, opts.init_config);
    gwm::dycore::NullTracerBoundaryUpdater null_tracer_boundary_updater;
    gwm::dycore::TracerBoundaryUpdater& tracer_boundary_updater =
        opts.enable_boundaries
            ? static_cast<gwm::dycore::TracerBoundaryUpdater&>(
                  prepared_tracer_boundary_updater)
            : static_cast<gwm::dycore::TracerBoundaryUpdater&>(
                  null_tracer_boundary_updater);
    gwm::dycore::LocalSplitExplicitFastMode local_fast_modes;
    NullFastModeIntegrator null_fast_modes;
    gwm::dycore::FastModeIntegrator& fast_modes =
        opts.enable_fast_modes
            ? static_cast<gwm::dycore::FastModeIntegrator&>(local_fast_modes)
            : static_cast<gwm::dycore::FastModeIntegrator&>(null_fast_modes);
    auto warm_rain_config = opts.warm_rain_config;
    warm_rain_config.dt = opts.step_config.dt;

    gwm::real sim_time = 0.0f;
    for (int step = 0; step < opts.steps; ++step) {
      trace_driver("step_" + std::to_string(step) + "_begin");
      if (opts.enable_boundaries) {
        prepared_boundary_updater.set_step_start_time(sim_time);
        prepared_tracer_boundary_updater.set_step_start_time(sim_time);
        trace_driver("step_" + std::to_string(step) + "_boundaries_configured");
      } else {
        trace_driver("step_" + std::to_string(step) + "_boundaries_disabled");
      }
      gwm::dycore::advance_dry_state_ssprk3(states, layout, metrics,
                                            opts.step_config, boundary_updater,
                                            fast_modes);
      trace_driver("step_" + std::to_string(step) + "_dry_done");
      maybe_write_phase_summary(step, sim_time, "dry", states, tracers, layout,
                                metrics, surface_precip_accum);
      if (opts.enable_tracer_transport) {
        gwm::dycore::advance_passive_tracers_ssprk3(
            tracers, states, layout, metrics, opts.step_config,
            tracer_boundary_updater);
        trace_driver("step_" + std::to_string(step) + "_tracers_done");
        maybe_write_phase_summary(step, sim_time, "tracers", states, tracers,
                                  layout, metrics, surface_precip_accum);
      } else {
        trace_driver("step_" + std::to_string(step) + "_tracers_skipped");
      }
      if (opts.enable_warm_rain) {
        gwm::physics::apply_warm_rain_microphysics(
            states, tracers, layout, metrics, warm_rain_config,
            &surface_precip_accum);
        trace_driver("step_" + std::to_string(step) + "_warm_rain_done");
        maybe_write_phase_summary(step, sim_time, "warm_rain", states, tracers,
                                  layout, metrics, surface_precip_accum);
      } else {
        trace_driver("step_" + std::to_string(step) + "_warm_rain_skipped");
      }
      sim_time += opts.step_config.dt;
    }

    GWM_CUDA_CHECK(cudaDeviceSynchronize());
    trace_driver("post_step_sync_done");

    const auto final =
        gwm::core::summarize_runtime_state(states, tracers, layout, metrics,
                                          &surface_precip_accum);
    trace_driver("final_summary_done");
    const auto summary_json =
        prepared_summary_json(runtime_case, opts.steps, opts.step_config.dt,
                              opts.step_config.fast_substeps,
                              opts.enable_boundaries,
                              opts.enable_fast_modes,
                              opts.enable_tracer_transport,
                              opts.enable_warm_rain,
                              opts.step_config.fast_update_horizontal_momentum,
                              opts.step_config.fast_update_density,
                              warm_rain_config,
                              startup_balance,
                              surface_runtime, initial, final);

    std::cout << summary_json << "\n";

    if (!opts.summary_json_path.empty()) {
      std::ofstream out(opts.summary_json_path, std::ios::binary);
      if (!out) {
        throw std::runtime_error("Failed to open summary path: " +
                                 opts.summary_json_path);
      }
      out << summary_json << "\n";
    }

    if (!opts.plan_view_json_path.empty()) {
      const auto plan_view = gwm::io::extract_runtime_plan_view(
          states, tracers, layout, metrics, "prepared_case", opts.steps,
          opts.step_config.dt, opts.plan_view_level, &surface_precip_accum);
      gwm::io::write_plan_view_bundle_json(plan_view, opts.plan_view_json_path);
    }

    GWM_CUDA_CHECK(cudaDeviceSynchronize());
    return EXIT_SUCCESS;
  } catch (const std::exception& ex) {
    std::cerr << "gwm_prepared_case_driver error: " << ex.what() << "\n";
    return EXIT_FAILURE;
  }
}
