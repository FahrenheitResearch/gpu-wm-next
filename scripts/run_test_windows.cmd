@echo off
setlocal EnableDelayedExpansion
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat" -arch=x64 -host_arch=x64 >nul
set CUDA_LAUNCH_BLOCKING=1
set BIN=%~1
shift
set ARGS=
:collect
if "%~1"=="" goto run
set ARGS=!ARGS! "%~1"
shift
goto collect
:run
call "%BIN%" !ARGS!
exit /b %ERRORLEVEL%
