param(
    [string]$BuildDir = "build-ninja",
    [string]$Generator = "Ninja",
    [string]$Configuration = "Release"
)

$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
$vsDevCmd = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat"
$cudaNvcc = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.0\bin\nvcc.exe"

if (-not (Test-Path $vsDevCmd)) {
    throw "VsDevCmd.bat not found at '$vsDevCmd'."
}

if (-not (Test-Path $cudaNvcc)) {
    throw "nvcc not found at '$cudaNvcc'."
}

$configure = @"
call "$vsDevCmd" -arch=x64 -host_arch=x64
cmake -S "$repoRoot" -B "$repoRoot\$BuildDir" -G "$Generator" -DCMAKE_BUILD_TYPE=$Configuration -DCMAKE_CUDA_COMPILER="$cudaNvcc"
if errorlevel 1 exit /b 1
cmake --build "$repoRoot\$BuildDir" --parallel
if errorlevel 1 exit /b 1
ctest --test-dir "$repoRoot\$BuildDir" --output-on-failure -C $Configuration
"@

cmd /c $configure
