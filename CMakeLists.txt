cmake_minimum_required(VERSION 3.27)

project(
  gpu_wm_next
  VERSION 0.1.0
  DESCRIPTION "Fresh-start GPU-first mesoscale weather model scaffold"
  LANGUAGES CXX CUDA
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(GWM_ENABLE_MPI "Enable MPI hooks when MPI is available" ON)

if(NOT CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES native)
endif()

include(CTest)
enable_testing()

find_package(CUDAToolkit REQUIRED)
find_package(MPI QUIET COMPONENTS CXX)

add_library(
  gwm_core
  STATIC
  src/comm/cartesian_topology.cpp
  src/comm/collectives.cpp
  src/comm/virtual_rank_layout.cpp
  src/comm/halo_exchange.cpp
  src/comm/mpi_runtime.cpp
  src/core/runtime_summary.cpp
  src/domain/idealized_domain_builder.cpp
  src/dycore/dry_core.cu
  src/dycore/boundary_conditions.cpp
  src/dycore/dry_diagnostics.cpp
  src/dycore/dry_fast_modes.cu
  src/dycore/dry_momentum_flux.cu
  src/dycore/dry_pressure_gradient.cu
  src/dycore/idealized_cases.cpp
  src/dycore/dry_transport.cu
  src/ingest/json_value.cpp
  src/ingest/prepared_case_init.cpp
  src/ingest/runtime_case.cpp
  src/ingest/source_catalog.cpp
  src/io/plan_view_output.cpp
  src/obsops/screen_obsops.cpp
  src/physics/warm_rain.cu
  src/surface/surface_layer_closure.cpp
  src/surface/surface_layer_exchange.cpp
)

target_include_directories(gwm_core PUBLIC include third_party)
target_compile_features(gwm_core PUBLIC cxx_std_20 cuda_std_20)
target_link_libraries(gwm_core PUBLIC CUDA::cudart)
set_target_properties(gwm_core PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

if(GWM_ENABLE_MPI AND MPI_CXX_FOUND)
  target_link_libraries(gwm_core PUBLIC MPI::MPI_CXX)
  target_compile_definitions(gwm_core PUBLIC GWM_HAVE_MPI=1)
else()
  target_compile_definitions(gwm_core PUBLIC GWM_HAVE_MPI=0)
endif()

if(GWM_ENABLE_MPI AND MPI_CXX_FOUND AND MPIEXEC_EXECUTABLE
   AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GwmAddMpiTest.cmake")
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/GwmAddMpiTest.cmake)
endif()

function(gwm_add_test target)
  add_executable(${target} ${ARGN})
  target_link_libraries(${target} PRIVATE gwm_core CUDA::cudart)
  target_include_directories(${target} PRIVATE tests)
  if(WIN32)
    add_test(
      NAME ${target}
      COMMAND
        ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_test_windows.cmd
        $<TARGET_FILE:${target}>
    )
  else()
    add_test(NAME ${target} COMMAND ${target})
  endif()
endfunction()

add_executable(gwm_idealized_driver src/core/gwm_idealized_driver.cpp)
target_link_libraries(gwm_idealized_driver PRIVATE gwm_core CUDA::cudart)
target_compile_features(gwm_idealized_driver PRIVATE cxx_std_20 cuda_std_20)

add_executable(gwm_prepared_case_driver src/core/gwm_prepared_case_driver.cpp)
target_link_libraries(gwm_prepared_case_driver PRIVATE gwm_core CUDA::cudart)
target_compile_features(gwm_prepared_case_driver PRIVATE cxx_std_20 cuda_std_20)

gwm_add_test(test_tracer_registry tests/unit/test_tracer_registry.cpp)
gwm_add_test(test_halo_exchange tests/unit/test_halo_exchange.cpp)
gwm_add_test(test_face_halo_exchange tests/unit/test_face_halo_exchange.cpp)
gwm_add_test(test_scalar_halo_buffers tests/unit/test_scalar_halo_buffers.cpp)
gwm_add_test(test_boundary_conditions tests/unit/test_boundary_conditions.cpp)
gwm_add_test(test_mpi_runtime tests/unit/test_mpi_runtime.cpp)
gwm_add_test(test_source_catalog tests/unit/test_source_catalog.cpp)
gwm_add_test(test_ingest_contracts tests/unit/test_ingest_contracts.cpp)
gwm_add_test(test_runtime_case tests/unit/test_runtime_case.cpp)
gwm_add_test(test_runtime_summary tests/unit/test_runtime_summary.cpp)
gwm_add_test(test_prepared_case_init tests/unit/test_prepared_case_init.cpp)
gwm_add_test(test_plan_view_output tests/unit/test_plan_view_output.cpp)
gwm_add_test(
  test_warm_rain_microphysics
  tests/unit/test_warm_rain_microphysics.cpp
)
gwm_add_test(test_dry_diagnostics tests/unit/test_dry_diagnostics.cpp)
gwm_add_test(test_dry_fast_modes tests/unit/test_dry_fast_modes.cpp)
gwm_add_test(test_dry_momentum_flux tests/unit/test_dry_momentum_flux.cpp)
gwm_add_test(test_dry_pressure_gradient tests/unit/test_dry_pressure_gradient.cpp)
gwm_add_test(
  test_slow_tendency_composition
  tests/unit/test_slow_tendency_composition.cpp
)
gwm_add_test(
  test_hybrid_height_metrics
  tests/unit/test_hybrid_height_metrics.cpp
)
gwm_add_test(
  test_idealized_domain_builder
  tests/unit/test_idealized_domain_builder.cpp
)
gwm_add_test(test_idealized_topography tests/unit/test_idealized_topography.cpp)
gwm_add_test(test_idealized_cases tests/unit/test_idealized_cases.cpp)
gwm_add_test(test_surface_state tests/unit/test_surface_state.cpp)
gwm_add_test(
  test_surface_static_properties
  tests/unit/test_surface_static_properties.cpp
)
gwm_add_test(test_surface_obsops tests/unit/test_surface_obsops.cpp)
gwm_add_test(
  test_surface_layer_closure
  tests/unit/test_surface_layer_closure.cpp
)
gwm_add_test(
  test_surface_tile_permutation
  tests/property/test_surface_tile_permutation.cpp
)
gwm_add_test(
  test_virtual_rank_equivalence
  tests/property/test_virtual_rank_equivalence.cpp
)
gwm_add_test(
  test_dry_virtual_rank_equivalence
  tests/property/test_dry_virtual_rank_equivalence.cpp
)
gwm_add_test(
  test_dry_momentum_virtual_rank_equivalence
  tests/property/test_dry_momentum_virtual_rank_equivalence.cpp
)
gwm_add_test(
  test_dry_fastmode_virtual_rank_equivalence
  tests/property/test_dry_fastmode_virtual_rank_equivalence.cpp
)
gwm_add_test(
  test_terrain_virtual_rank_equivalence
  tests/property/test_terrain_virtual_rank_equivalence.cpp
)

if(GWM_ENABLE_MPI AND MPI_CXX_FOUND AND MPIEXEC_EXECUTABLE
   AND COMMAND gwm_add_mpi_test)
  gwm_add_mpi_test(
    test_mpi_cartesian_context
    4
    tests/unit/test_mpi_cartesian_context.cpp
  )
  gwm_add_mpi_test(
    test_mpi_scalar_halo_exchange
    4
    tests/unit/test_mpi_scalar_halo_exchange.cpp
  )
  gwm_add_mpi_test(
    test_mpi_face_halo_exchange
    4
    tests/unit/test_mpi_face_halo_exchange.cpp
  )
  gwm_add_mpi_test(
    test_mpi_face_interface_sync
    4
    tests/unit/test_mpi_face_interface_sync.cpp
  )
  gwm_add_mpi_test(
    test_mpi_halo_mismatch_rejects
    4
    tests/unit/test_mpi_halo_mismatch_rejects.cpp
  )
  gwm_add_mpi_test(
    test_mpi_halo_virtual_rank_equivalence
    4
    tests/property/test_mpi_halo_virtual_rank_equivalence.cpp
  )
endif()
gwm_add_test(test_constant_state tests/regression/test_constant_state.cpp)
gwm_add_test(
  test_dry_constant_state_bundle
  tests/regression/test_dry_constant_state_bundle.cpp
)
gwm_add_test(
  test_hydrostatic_rest
  tests/regression/test_hydrostatic_rest.cpp
)
gwm_add_test(
  test_terrain_hydrostatic_rest
  tests/regression/test_terrain_hydrostatic_rest.cpp
)
gwm_add_test(
  test_buoyancy_response
  tests/regression/test_buoyancy_response.cpp
)
gwm_add_test(
  test_horizontal_pressure_response
  tests/regression/test_horizontal_pressure_response.cpp
)
gwm_add_test(
  test_density_current_evolution
  tests/regression/test_density_current_evolution.cpp
)
gwm_add_test(
  test_acoustic_pulse
  tests/regression/test_acoustic_pulse.cpp
)
gwm_add_test(
  test_mountain_wave_evolution
  tests/regression/test_mountain_wave_evolution.cpp
)
gwm_add_test(
  test_surface_layer_neutral_reference
  tests/regression/test_surface_layer_neutral_reference.cpp
)
gwm_add_test(
  test_scalar_mass_conservation
  tests/regression/test_scalar_mass_conservation.cpp
)
gwm_add_test(
  test_warm_rain_closed_box
  tests/regression/test_warm_rain_closed_box.cpp
)
gwm_add_test(
  test_warm_rain_summary_closure
  tests/regression/test_warm_rain_summary_closure.cpp
)
