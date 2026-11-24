[CmdletBinding()]
param(
    [ValidateSet('native', 'mkl', 'mkl-mumps')]
    [string[]]$Profile = @('native'),
    [ValidateSet('Debug', 'Release')]
    [string[]]$Configuration = @('Release'),
    [string]$BuildRoot,
    [string]$DependencyCache,
    [string]$MklDir,
    [string]$NinjaPath,
    [string]$SmokeScript,
    [ValidateRange(1, 256)]
    [int]$Jobs = 4
)

$ErrorActionPreference = 'Stop'
$source = (Get-Item -LiteralPath (Join-Path $PSScriptRoot '..')).FullName
if (-not $BuildRoot) {
    $BuildRoot = Join-Path $source 'build/local-regression'
}
$BuildRoot = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($BuildRoot)

$dependencyArgs = @()
if ($DependencyCache) {
    $dependencyArgs += '-C', (Resolve-Path -LiteralPath $DependencyCache).ProviderPath
}
if (-not $MklDir -and $env:MKLROOT) {
    $candidate = Join-Path $env:MKLROOT 'lib/cmake/mkl'
    if (Test-Path -LiteralPath (Join-Path $candidate 'MKLConfig.cmake')) {
        $MklDir = $candidate
    }
}
if (-not $NinjaPath -and $env:VSINSTALLDIR -and -not (Get-Command ninja -ErrorAction SilentlyContinue)) {
    $candidate = Join-Path $env:VSINSTALLDIR 'Common7/IDE/CommonExtensions/Microsoft/CMake/Ninja/ninja.exe'
    if (Test-Path -LiteralPath $candidate) {
        $NinjaPath = $candidate
    }
}
if ($NinjaPath) {
    $dependencyArgs += "-DCMAKE_MAKE_PROGRAM:FILEPATH=$((Resolve-Path -LiteralPath $NinjaPath).ProviderPath)"
}
if ($SmokeScript) {
    $SmokeScript = (Resolve-Path -LiteralPath $SmokeScript).ProviderPath
}

foreach ($backend in $Profile) {
    $useMkl = if ($backend -eq 'native') { 'OFF' } else { 'ON' }
    $useMumps = if ($backend -eq 'mkl-mumps') { 'ON' } else { 'OFF' }
    $mklArgs = @()
    if ($useMkl -eq 'ON' -and $MklDir) {
        $mklArgs += "-DMKL_DIR:PATH=$((Resolve-Path -LiteralPath $MklDir).ProviderPath)"
    }
    foreach ($config in $Configuration) {
        $build = Join-Path $BuildRoot "$backend-$config"
        Write-Host "Building $backend / $config in $build"
        $profileArgs = @(
            "-DCMAKE_BUILD_TYPE=$config"
            "-DBPFEM_USE_MKL=$useMkl"
            "-DBPFEM_REQUIRE_MKL=$useMkl"
            "-DBPFEM_USE_MUMPS=$useMumps"
        )
        & cmake -S $source -B $build -G Ninja @dependencyArgs @mklArgs @profileArgs
        if ($LASTEXITCODE -ne 0) { throw "CMake configure failed for $backend / $config" }

        & cmake --build $build --parallel $Jobs
        if ($LASTEXITCODE -ne 0) { throw "Build failed for $backend / $config" }

        & ctest --test-dir $build --output-on-failure --timeout 120 --no-tests=error
        if ($LASTEXITCODE -ne 0) { throw "CTest failed for $backend / $config" }

        if ($SmokeScript) {
            $exeName = if ($env:OS -eq 'Windows_NT') { 'bp_fem_solver.exe' } else { 'bp_fem_solver' }
            $solver = Join-Path $build $exeName
            $smokeOutput = Join-Path $build 'smoke'
            New-Item -ItemType Directory -Path $smokeOutput -Force | Out-Null
            $LASTEXITCODE = 0
            & $SmokeScript -SolverPath $solver -OutputDir $smokeOutput -BackendProfile $backend
            if (-not $? -or $LASTEXITCODE -ne 0) { throw "CLI smoke test failed for $backend / $config" }
        }
    }
}
