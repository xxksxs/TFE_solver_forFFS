#pragma once

#include <string>

namespace fem {

class Logger;
struct Options;

// Snapshot of build / runtime / host environment, captured once per run and
// emitted to the log header. The contents are descriptive only; nothing in
// here changes solver behavior.
//
// Fields that cannot be queried on a given platform fall back to "n/a";
// the caller is expected to print every field regardless so the log layout
// stays uniform across machines.
struct RunEnvironment {
    std::string solverVersion;       // e.g. "bp_fem_solver 2026R1"
    std::string buildString;         // compiler + config + feature defines
    std::string runtimeString;       // MKL / threads info if available
    std::string hostName;            // computer name
    std::string osDescription;       // OS family + version
    std::string cpuDescription;      // CPU brand string + logical core count
    std::string memoryDescription;   // total / available physical memory
    std::string commandLine;         // argv joined with spaces
    std::string workingDirectory;    // CWD at process start
};

// Capture all fields. Safe to call once at start-of-run; cheap.
RunEnvironment captureRunEnvironment(int argc, char** argv);

// Emit `[env] key  value` lines for every captured field. Intended to be
// called immediately after Logger::banner() and Logger::attachLogFile(),
// before the first phase().
void writeEnvHeader(Logger& log, const RunEnvironment& env);

}  // namespace fem
