#include <cstdlib>
#include <fstream>
#include <iostream>
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
  int steps = 1;
  std::string summary_json_path;
  std::string plan_view_json_path;
  int plan_view_level = 0;
};

bool driver_trace_enabled() {
  const char* value = std::getenv("GWM_TRACE_DRIVER");
  if (value == nullptr) {
    return false;
  }
  return std::string(value) != "0" && std::string(value) != "false";
}

void trace_driver(const std::string& message) {
  if (!driver_trace_enabled()) {
    return;
  }
  std::cerr << "[gwm_prepared_case_driver] " << message << std::endl;
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
    } else if (arg == "--dt") {
      opts.step_config.dt = parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--fast-substeps") {
      opts.step_config.fast_substeps =
          parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--gravity") {
      opts.step_config.gravity = parse_value<float>(require_value(i, argc, argv));
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
    gwm::real dt, int fast_substeps,
    const gwm::surface::SurfaceRuntimeInitResult& surface_runtime,
    const gwm::core::RuntimeStateSummary& initial,
    const gwm::core::RuntimeStateSummary& final) {
  const double elapsed_seconds =
      static_cast<double>(steps) * static_cast<double>(dt);
  const double elapsed_hours = elapsed_seconds / 3600.0;
  std::ostringstream oss;
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
        runtime_case.analysis, layout, "prepared_case");
    trace_driver("dry_states_initialized");
    auto tracers = gwm::ingest::make_warm_rain_tracers_from_analysis(
        runtime_case.analysis, states, layout, "prepared_case_tracer");
    trace_driver("warm_rain_tracers_initialized");

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

    gwm::ingest::PreparedCaseBoundaryUpdater boundary_updater(
        runtime_case.analysis, runtime_case.boundary_cache, opts.init_config);
    gwm::ingest::PreparedCaseTracerBoundaryUpdater tracer_boundary_updater(
        runtime_case.analysis, runtime_case.boundary_cache, opts.init_config);
    gwm::dycore::LocalSplitExplicitFastMode fast_modes;
    gwm::physics::WarmRainConfig warm_rain_config{};
    warm_rain_config.dt = opts.step_config.dt;
    warm_rain_config.rain_terminal_velocity = 12.0f;

    gwm::real sim_time = 0.0f;
    for (int step = 0; step < opts.steps; ++step) {
      trace_driver("step_" + std::to_string(step) + "_begin");
      boundary_updater.set_step_start_time(sim_time);
      tracer_boundary_updater.set_step_start_time(sim_time);
      trace_driver("step_" + std::to_string(step) + "_boundaries_configured");
      gwm::dycore::advance_dry_state_ssprk3(states, layout, metrics,
                                            opts.step_config, boundary_updater,
                                            fast_modes);
      trace_driver("step_" + std::to_string(step) + "_dry_done");
      gwm::dycore::advance_passive_tracers_ssprk3(
          tracers, states, layout, metrics, opts.step_config,
          tracer_boundary_updater);
      trace_driver("step_" + std::to_string(step) + "_tracers_done");
      gwm::physics::apply_warm_rain_microphysics(
          states, tracers, layout, metrics, warm_rain_config,
          &surface_precip_accum);
      trace_driver("step_" + std::to_string(step) + "_warm_rain_done");
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
                              opts.step_config.fast_substeps, surface_runtime,
                              initial, final);

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
