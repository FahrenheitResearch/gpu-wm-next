#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gwm/domain/idealized_domain_builder.hpp"
#include "gwm/dycore/dry_core.hpp"
#include "gwm/dycore/dry_diagnostics.hpp"
#include "gwm/dycore/idealized_cases.hpp"
#include "gwm/io/plan_view_output.hpp"

namespace {

using gwm::domain::RectilinearDomainConfig;
using gwm::dycore::MountainWaveBackgroundConfig;
using gwm::dycore::ThermoBubbleConfig;

struct DriverOptions {
  std::string case_kind = "density_current";
  RectilinearDomainConfig domain{};
  ThermoBubbleConfig bubble{};
  MountainWaveBackgroundConfig mountain{};
  gwm::dycore::DryStepperConfig step_config{};
  int steps = 10;
  std::string summary_json_path;
  std::string plan_view_json_path;
  int plan_view_level = 0;
};

bool parse_bool(const std::string& value) {
  if (value == "true" || value == "1") {
    return true;
  }
  if (value == "false" || value == "0") {
    return false;
  }
  throw std::runtime_error("Invalid boolean value: " + value);
}

template <typename T>
T parse_value(const std::string& value);

template <>
int parse_value<int>(const std::string& value) {
  return std::stoi(value);
}

template <>
double parse_value<double>(const std::string& value) {
  return std::stod(value);
}

template <>
float parse_value<float>(const std::string& value) {
  return std::stof(value);
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
    if (arg == "--case") {
      opts.case_kind = require_value(i, argc, argv);
    } else if (arg == "--nx") {
      opts.domain.nx = parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--ny") {
      opts.domain.ny = parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--nz") {
      opts.domain.nz = parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--halo") {
      opts.domain.halo = parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--ranks-x") {
      opts.domain.ranks_x = parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--ranks-y") {
      opts.domain.ranks_y = parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--dx") {
      opts.domain.dx = parse_value<double>(require_value(i, argc, argv));
    } else if (arg == "--dy") {
      opts.domain.dy = parse_value<double>(require_value(i, argc, argv));
    } else if (arg == "--z-top") {
      opts.domain.z_top = parse_value<double>(require_value(i, argc, argv));
    } else if (arg == "--periodic-x") {
      opts.domain.periodic_x = parse_bool(require_value(i, argc, argv));
    } else if (arg == "--periodic-y") {
      opts.domain.periodic_y = parse_bool(require_value(i, argc, argv));
    } else if (arg == "--terrain-kind") {
      const auto value = require_value(i, argc, argv);
      if (value == "flat") {
        opts.domain.terrain_kind = gwm::domain::TerrainProfileKind::Flat;
      } else if (value == "cosine_mountain") {
        opts.domain.terrain_kind =
            gwm::domain::TerrainProfileKind::CosineMountain;
      } else {
        throw std::runtime_error("Unsupported terrain kind: " + value);
      }
    } else if (arg == "--terrain-center-x") {
      opts.domain.mountain_center_x =
          parse_value<double>(require_value(i, argc, argv));
    } else if (arg == "--terrain-center-y") {
      opts.domain.mountain_center_y =
          parse_value<double>(require_value(i, argc, argv));
    } else if (arg == "--terrain-half-width-x") {
      opts.domain.mountain_half_width_x =
          parse_value<double>(require_value(i, argc, argv));
    } else if (arg == "--terrain-half-width-y") {
      opts.domain.mountain_half_width_y =
          parse_value<double>(require_value(i, argc, argv));
    } else if (arg == "--terrain-height") {
      opts.domain.mountain_height =
          parse_value<double>(require_value(i, argc, argv));
    } else if (arg == "--rho-background") {
      opts.bubble.rho_background =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--theta-background") {
      opts.bubble.theta_background =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--u-background") {
      const auto value = parse_value<float>(require_value(i, argc, argv));
      opts.bubble.u_background = value;
      opts.mountain.u_background = value;
    } else if (arg == "--v-background") {
      const auto value = parse_value<float>(require_value(i, argc, argv));
      opts.bubble.v_background = value;
      opts.mountain.v_background = value;
    } else if (arg == "--w-background") {
      const auto value = parse_value<float>(require_value(i, argc, argv));
      opts.bubble.w_background = value;
      opts.mountain.w_background = value;
    } else if (arg == "--theta-perturbation") {
      opts.bubble.theta_perturbation =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--center-x-fraction") {
      opts.bubble.center_x_fraction =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--center-y-fraction") {
      opts.bubble.center_y_fraction =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--center-z-fraction") {
      opts.bubble.center_z_fraction =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--radius-x-fraction") {
      opts.bubble.radius_x_fraction =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--radius-y-fraction") {
      opts.bubble.radius_y_fraction =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--radius-z-fraction") {
      opts.bubble.radius_z_fraction =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--rho-surface") {
      opts.mountain.rho_surface =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--theta-ref") {
      opts.mountain.theta_ref = parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--density-scale-height") {
      opts.mountain.density_scale_height =
          parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--steps") {
      opts.steps = parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--dt") {
      opts.step_config.dt = parse_value<float>(require_value(i, argc, argv));
    } else if (arg == "--fast-substeps") {
      opts.step_config.fast_substeps =
          parse_value<int>(require_value(i, argc, argv));
    } else if (arg == "--gravity") {
      opts.step_config.gravity =
          parse_value<float>(require_value(i, argc, argv));
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
}

std::vector<gwm::dycore::DryState> make_case_state(
    const DriverOptions& opts, const gwm::domain::IdealizedDomain& domain) {
  if (opts.case_kind == "warm_bubble") {
    return gwm::dycore::make_warm_bubble_state(domain.layout, domain.metrics,
                                               opts.bubble, "warm_bubble");
  }
  if (opts.case_kind == "density_current") {
    return gwm::dycore::make_density_current_state(domain.layout, domain.metrics,
                                                   opts.bubble,
                                                   "density_current");
  }
  if (opts.case_kind == "hydrostatic_rest") {
    return gwm::dycore::make_hydrostatic_rest_state(
        domain.layout, domain.metrics, opts.mountain.rho_surface,
        opts.mountain.theta_ref, opts.mountain.density_scale_height,
        "hydrostatic_rest");
  }
  if (opts.case_kind == "mountain_wave_background") {
    return gwm::dycore::make_mountain_wave_background_state(
        domain.layout, domain.metrics, domain, opts.mountain,
        "mountain_wave");
  }
  throw std::runtime_error("Unsupported case kind: " + opts.case_kind);
}

std::string driver_summary_json(const std::string& case_kind, int steps,
                                gwm::real dt,
                                const gwm::dycore::DryStateSummary& initial,
                                const gwm::dycore::DryStateSummary& final) {
  std::ostringstream oss;
  oss << "{\n";
  oss << "  \"case\": \"" << case_kind << "\",\n";
  oss << "  \"steps\": " << steps << ",\n";
  oss << "  \"dt\": " << dt << ",\n";
  oss << "  \"initial\": "
      << gwm::dycore::summary_to_json(initial, "    ") << ",\n";
  oss << "  \"final\": " << gwm::dycore::summary_to_json(final, "    ") << "\n";
  oss << "}";
  return oss.str();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    DriverOptions opts{};
    parse_args(argc, argv, opts);

    const auto domain = gwm::domain::build_rectilinear_domain(opts.domain);
    auto states = make_case_state(opts, domain);
    const auto initial = gwm::dycore::summarize_dry_states(states);

    gwm::dycore::NullBoundaryUpdater boundary;
    gwm::dycore::LocalSplitExplicitFastMode fast_modes;
    for (int step = 0; step < opts.steps; ++step) {
      gwm::dycore::advance_dry_state_ssprk3(states, domain.layout,
                                            domain.metrics, opts.step_config,
                                            boundary, fast_modes);
    }
    const auto final = gwm::dycore::summarize_dry_states(states);

    const auto json =
        driver_summary_json(opts.case_kind, opts.steps, opts.step_config.dt,
                            initial, final);
    std::cout << json << "\n";

    if (!opts.summary_json_path.empty()) {
      std::ofstream out(opts.summary_json_path, std::ios::binary);
      if (!out) {
        throw std::runtime_error("Failed to open summary path: " +
                                 opts.summary_json_path);
      }
      out << json << "\n";
    }

    if (!opts.plan_view_json_path.empty()) {
      const auto plan_view = gwm::io::extract_dry_plan_view(
          states, domain.layout, domain.metrics, opts.case_kind, opts.steps,
          opts.step_config.dt, opts.plan_view_level);
      gwm::io::write_plan_view_bundle_json(plan_view, opts.plan_view_json_path);
    }
    return EXIT_SUCCESS;
  } catch (const std::exception& ex) {
    std::cerr << "gwm_idealized_driver error: " << ex.what() << "\n";
    return EXIT_FAILURE;
  }
}
