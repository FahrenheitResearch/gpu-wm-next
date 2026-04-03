@echo off
setlocal

call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat" -arch=x64 -host_arch=x64 >nul

set ROOT=%~dp0..
set BUILD=%ROOT%\build-ninja

"%BUILD%\test_idealized_domain_builder.exe"
if errorlevel 1 exit /b 1

"%BUILD%\test_dry_virtual_rank_equivalence.exe"
if errorlevel 1 exit /b 1

"%BUILD%\test_dry_momentum_virtual_rank_equivalence.exe"
if errorlevel 1 exit /b 1

"%BUILD%\test_dry_momentum_flux.exe"
if errorlevel 1 exit /b 1

"%BUILD%\test_dry_fastmode_virtual_rank_equivalence.exe"
if errorlevel 1 exit /b 1

"%BUILD%\test_dry_fast_modes.exe"
if errorlevel 1 exit /b 1

"%BUILD%\test_dry_constant_state_bundle.exe"
if errorlevel 1 exit /b 1

"%BUILD%\test_hydrostatic_rest.exe"
if errorlevel 1 exit /b 1

"%BUILD%\test_buoyancy_response.exe"
if errorlevel 1 exit /b 1

"%BUILD%\test_horizontal_pressure_response.exe"
if errorlevel 1 exit /b 1

"%BUILD%\test_density_current_evolution.exe"
if errorlevel 1 exit /b 1

"%BUILD%\test_acoustic_pulse.exe"
if errorlevel 1 exit /b 1

exit /b 0
