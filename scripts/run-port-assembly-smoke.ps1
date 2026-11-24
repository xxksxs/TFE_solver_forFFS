[CmdletBinding()]
param(
    [Parameter(Mandatory)][string]$SolverPath,
    [Parameter(Mandatory)][string]$OutputDir,
    [Parameter(Mandatory)][ValidateSet('native', 'mkl', 'mkl-mumps')][string]$BackendProfile,
    [string]$BenchmarkRoot = (Join-Path $PSScriptRoot '../../../benchmark')
)

$ErrorActionPreference = 'Stop'
if ($BackendProfile -eq 'native') {
    throw 'This regression uses PARDISO for its Direct reference. Select -Profile mkl or mkl-mumps.'
}

$solver = (Resolve-Path -LiteralPath $SolverPath).ProviderPath
$exeSuffix = [IO.Path]::GetExtension($solver)
$probe = Join-Path (Split-Path -Parent $solver) ("io_parser_tests" + $exeSuffix)
$benchmark = (Resolve-Path -LiteralPath $BenchmarkRoot).ProviderPath
$node = (Get-Command node -CommandType Application -ErrorAction Stop | Select-Object -First 1).Source
# The smoke script requires a new output directory for each run.
# Use a unique child of the build's smoke directory to preserve earlier evidence.
$run = Join-Path $OutputDir ("port-assembly-" + [Guid]::NewGuid().ToString('N'))
$nodeArgs = @(
    (Join-Path $PSScriptRoot 'port_assembly_smoke.mjs')
    '--solver', $solver, '--probe', $probe, '--benchmark', $benchmark, '--out', $run
)
& $node @nodeArgs
if ($LASTEXITCODE -ne 0) {
    throw "Port assembly smoke failed (exit code $LASTEXITCODE). Evidence: $run"
}
