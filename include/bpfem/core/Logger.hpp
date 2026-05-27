#pragma once

#include <chrono>
#include <deque>
#include <fstream>
#include <iosfwd>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace fem {

// Per-run process memory snapshot (resident set / working set, peak working
// set). Both fields are bytes; sentinel value 0 means "not available on this
// platform / could not read" and is treated as a missing measurement.
struct MemoryStats {
    std::size_t currentBytes = 0;
    std::size_t peakBytes = 0;
};

MemoryStats queryProcessMemory();

// Format a byte count using kB / MB / GB units with 2 decimals. "1.50 GB",
// "234.50 MB", "12.30 kB".
std::string formatBytes(std::size_t bytes);

// Tee logger: every banner / info / phase line is written both to a console
// stream (typically std::cout) and -- when configured -- to a log file in the
// run's output directory. Phase boundaries also accumulate timing and memory
// snapshots so a final summary() call can print a cumulative table.
//
// Adds (Phase 1++ extensions):
//   - structured tag-style lines: env(), stat(), warn(), error()
//   - context stack used by writeCrashDump() to identify what was running
//     at fatal time
//   - a recent-line ring buffer so writeCrashDump() can flush the last few
//     hundred lines into a separate run.crash.log
//
// Usage:
//   Logger log(std::cout);
//   log.attachLogFile(outDir / "run.log");
//   log.banner();
//   log.env("solver", "bp_fem_solver 2026R1");        // env header lines
//   log.phase("Project import");
//     ... log.info(...) ...
//   log.summary();                                     // table at end
class Logger {
public:
    explicit Logger(std::ostream& consoleStream);
    ~Logger();

    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    // Open `path` for truncate-and-write and tee subsequent output to it.
    // If the file cannot be opened, an info() line records the failure and
    // only the console destination remains.
    void attachLogFile(const std::string& path);

    // ---- Lines that go to console + log file ----
    void banner();
    void info(const std::string& message);
    void warn(const std::string& message);
    void error(const std::string& message);

    // [env] key=value lines for the run header. Use early in main, before
    // the first phase().
    void env(const std::string& key, const std::string& value);

    // [stat] tagged structured lines for tabular output (S parameter table,
    // sub-step timing, etc.). Caller pre-formats the row. The leading "[stat] "
    // prefix is added automatically.
    void stat(const std::string& row);

    // Begin a new named phase. Closes the previous phase, recording its
    // elapsed wall time + end-of-phase memory snapshot, then prints the
    // dashed header + initial memory snapshot of the new phase.
    void phase(const std::string& message);

    // Per-frequency / per-step context tags. Used by writeCrashDump() to
    // tell the user which point was being solved when the failure occurred.
    // The Application / Sweep code calls pushContext("freq_GHz", "41.5")
    // before each frequency, popContext() after.
    void pushContext(const std::string& key, const std::string& value);
    void popContext();

    // Print the per-phase + total summary table.
    void summary();

    // Dump the last bufferedLines + active phase + active context into a
    // separate .crash.log next to the main log. Safe to call inside
    // catch(std::exception). Returns true iff the file was successfully
    // written.
    bool writeCrashDump(const std::string& crashLogPath, const std::string& whatMessage);

private:
    struct PhaseRecord {
        std::string name;
        std::chrono::steady_clock::time_point start;
        double elapsedSec = 0.0;
        MemoryStats startMem{};
        MemoryStats endMem{};
        bool closed = false;
    };

    void writeLine(const std::string& line);
    void closeCurrentPhase();
    static std::string timestamp();

    std::ostream& console_;
    std::unique_ptr<std::ofstream> file_;
    std::vector<PhaseRecord> phases_;
    std::vector<std::pair<std::string, std::string>> context_;
    std::deque<std::string> recentLines_;  // ring buffer for crash dump
    static constexpr std::size_t kRecentCapacity = 512;
    std::chrono::steady_clock::time_point runStart_{};
    bool runStarted_ = false;
    bool summaryEmitted_ = false;
};

}  // namespace fem
