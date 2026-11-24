#pragma once

#include "bpfem/core/Types.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace fem {

struct RunEnvironment;
struct Options;

// Per-phase timing + memory record for the structured JSON sidecar. Mirrors
// what Logger::summary() prints, but in a machine-readable form so CI / the
// post-processor exe can parse runs without grepping text logs.
struct PhaseSummary {
    std::string name;
    double elapsedSec = 0.0;
    double endCurrentMb = 0.0;
    double endPeakMb = 0.0;
};

// Per-frequency-point stats: assemble + linear-solve sub-step timing in
// addition to S parameters. Populated by Application's sweep loop wrapper.
struct FrequencyStat {
    double frequencyHz = 0.0;
    double assembleSec = 0.0;
    double solveSec = 0.0;
    double residual = 0.0;
    double absoluteResidual = 0.0;
    std::string linearSolveStatus;
    int iterations = 0;
    double s11Db = 0.0;
    double s21Db = 0.0;
};

// Aggregate run report. Built incrementally over the run, then dumped to
// `<outDir>/run.json` near the end. Kept hand-written (no third-party JSON
// dep) -- the format is small enough.
struct RunReport {
    int schema = 2;
    std::string status = "fatal";
    std::string solverVersion;
    std::string buildString;
    std::string runtimeString;
    std::string hostName;
    std::string osDescription;
    std::string cpuDescription;
    std::string memoryDescription;
    std::string commandLine;
    std::string workingDirectory;
    std::string sweepStrategyName;
    std::string linearSolverBackend;
    std::string outputDirectory;
    std::string aedtSha256;
    std::string meshSha256;
    std::string executableSha256;
    int mklMaxThreads = 0;
    int mklPardisoMaxThreads = 0;
    std::string ompNumThreads;

    // Mirrors Options. Stringified to keep the JSON writer trivial; CI can
    // still grep / regex over these.
    std::vector<std::pair<std::string, std::string>> options;

    std::vector<PhaseSummary> phases;
    std::vector<FrequencyStat> sweep;

    double totalElapsedSec = 0.0;
    double peakMemoryMb = 0.0;

    double offlineBuildSec = 0.0;
    double portLinearizationSec = 0.0;
    double lanczosOperatorSec = 0.0;
    double poleDecompositionSec = 0.0;
    double orthogonalizationSec = 0.0;
    double romProjectionSec = 0.0;
    double onlineSweepSec = 0.0;

    std::uint64_t symbolicAnalysisCount = 0;
    std::uint64_t numericFactorizationCount = 0;
    std::uint64_t factorizedRhsSolveCount = 0;
    std::uint64_t factorizedSolveCallCount = 0;
    std::uint64_t batchRhsMax = 0;
    std::uint64_t estimatedSolverPeakKiB = 0;
    double symbolicAnalysisSec = 0.0;
    double numericFactorizationSec = 0.0;
    double factorizedRhsSolveSec = 0.0;
};

// Populate the report's environment / options sections from a captured
// RunEnvironment plus the parsed Options. Does not touch phases / sweep.
void fillReportFromEnvironment(RunReport& report, const RunEnvironment& env, const Options& options);
void captureRunInputs(RunReport& report, const Options& options, const char* executable);
void verifyRunInputs(const RunReport& report, const Options& options);

// Write the report to disk. Returns true on success. Hand-rolled JSON; no
// schema validation.
bool writeRunReport(const std::filesystem::path& path, const RunReport& report);

// Write a compact timing/memory report intended for algorithm comparisons.
// This duplicates only the metrics most useful for comparing direct / ALPS /
// AWE-family runs, while run.json remains the full run sidecar.
bool writeTimingReport(const std::filesystem::path& path, const RunReport& report);

}  // namespace fem
