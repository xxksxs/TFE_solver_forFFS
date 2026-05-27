#include "bpfem/core/Logger.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <ios>
#include <ostream>
#include <sstream>
#include <utility>

#if defined(_WIN32)
  #define WIN32_LEAN_AND_MEAN
  #define NOMINMAX
  #include <windows.h>
  #include <psapi.h>
#elif defined(__linux__)
  #include <unistd.h>
  #include <fstream>
#elif defined(__APPLE__)
  #include <mach/mach.h>
#endif

namespace fem {

namespace {

std::string padRight(const std::string& s, std::size_t w) {
    if (s.size() >= w) {
        return s;
    }
    return s + std::string(w - s.size(), ' ');
}

}  // namespace

MemoryStats queryProcessMemory() {
    MemoryStats out;
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS info;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info))) {
        out.currentBytes = static_cast<std::size_t>(info.WorkingSetSize);
        out.peakBytes = static_cast<std::size_t>(info.PeakWorkingSetSize);
    }
#elif defined(__linux__)
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        const auto parseKb = [&](std::size_t& target) {
            const auto colon = line.find(':');
            if (colon == std::string::npos) return;
            std::istringstream iss(line.substr(colon + 1));
            long long kb = 0;
            std::string unit;
            iss >> kb >> unit;
            if (kb > 0) {
                target = static_cast<std::size_t>(kb) * 1024U;
            }
        };
        if (line.rfind("VmRSS:", 0) == 0) {
            parseKb(out.currentBytes);
        } else if (line.rfind("VmHWM:", 0) == 0) {
            parseKb(out.peakBytes);
        }
    }
#elif defined(__APPLE__)
    mach_task_basic_info info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                  reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS) {
        out.currentBytes = static_cast<std::size_t>(info.resident_size);
        out.peakBytes = out.currentBytes;
    }
#endif
    return out;
}

std::string formatBytes(std::size_t bytes) {
    constexpr double kKilo = 1024.0;
    constexpr double kMega = 1024.0 * 1024.0;
    constexpr double kGiga = 1024.0 * 1024.0 * 1024.0;
    std::ostringstream out;
    out << std::fixed << std::setprecision(2);
    const double b = static_cast<double>(bytes);
    if (b >= kGiga) {
        out << (b / kGiga) << " GB";
    } else if (b >= kMega) {
        out << (b / kMega) << " MB";
    } else if (b >= kKilo) {
        out << (b / kKilo) << " kB";
    } else {
        out << bytes << " B";
    }
    return out.str();
}

Logger::Logger(std::ostream& consoleStream) : console_(consoleStream) {}

Logger::~Logger() {
    if (runStarted_ && !summaryEmitted_) {
        summary();
    }
}

void Logger::attachLogFile(const std::string& path) {
    file_ = std::make_unique<std::ofstream>(path, std::ios::out | std::ios::trunc);
    if (!file_->good()) {
        file_.reset();
        info("Could not open log file: " + path);
        return;
    }
    info("Tee logging to file: " + path);
}

void Logger::banner() {
    const std::string lines[] = {
        "",
        "============================================================",
        "  BP-FEM Frequency Domain Solver 2026R1",
        "  Driven-modal waveguide workflow | AEDT + NGMesh frontend",
        "============================================================",
    };
    for (const auto& l : lines) {
        writeLine(l);
    }
}

void Logger::info(const std::string& message) {
    writeLine(timestamp() + " [info] " + message);
}

void Logger::warn(const std::string& message) {
    writeLine(timestamp() + " [warn] " + message);
}

void Logger::error(const std::string& message) {
    writeLine(timestamp() + " [error] " + message);
}

void Logger::env(const std::string& key, const std::string& value) {
    // Fixed-width key column makes the env block scannable. Width 14 fits
    // every key we currently emit ("solver", "build", "runtime", "host",
    // "cpu", "memory", "cmdline", "cwd", ...) without wrapping.
    constexpr std::size_t kKeyWidth = 14;
    writeLine("[env] " + padRight(key, kKeyWidth) + value);
}

void Logger::stat(const std::string& row) {
    writeLine("[stat] " + row);
}

void Logger::phase(const std::string& message) {
    if (!runStarted_) {
        runStart_ = std::chrono::steady_clock::now();
        runStarted_ = true;
    }

    closeCurrentPhase();

    PhaseRecord rec;
    rec.name = message;
    rec.start = std::chrono::steady_clock::now();
    rec.startMem = queryProcessMemory();
    phases_.push_back(rec);

    writeLine("");
    writeLine("--- " + message + " ---");
    if (rec.startMem.currentBytes > 0) {
        std::ostringstream oss;
        oss << timestamp() << " [info] memory: current="
            << formatBytes(rec.startMem.currentBytes)
            << ", peak=" << formatBytes(rec.startMem.peakBytes);
        writeLine(oss.str());
    }
}

void Logger::pushContext(const std::string& key, const std::string& value) {
    context_.emplace_back(key, value);
}

void Logger::popContext() {
    if (!context_.empty()) {
        context_.pop_back();
    }
}

void Logger::summary() {
    if (summaryEmitted_) {
        return;
    }
    closeCurrentPhase();
    summaryEmitted_ = true;

    const auto totalEnd = std::chrono::steady_clock::now();
    const double totalSec = runStarted_
        ? std::chrono::duration<double>(totalEnd - runStart_).count()
        : 0.0;
    const auto finalMem = queryProcessMemory();

    writeLine("");
    writeLine("--- Run summary ---");

    constexpr std::size_t kPhaseW = 28;
    constexpr std::size_t kTimeW = 12;
    constexpr std::size_t kMemW = 14;
    std::ostringstream header;
    header << "  " << padRight("phase", kPhaseW)
           << padRight("elapsed", kTimeW)
           << padRight("end-current", kMemW)
           << padRight("end-peak", kMemW);
    writeLine(header.str());
    writeLine("  " + std::string(kPhaseW + kTimeW + 2 * kMemW - 2, '-'));

    for (const auto& p : phases_) {
        std::ostringstream row;
        row << "  " << padRight(p.name, kPhaseW);
        std::ostringstream tcol;
        tcol << std::fixed << std::setprecision(3) << p.elapsedSec << " s";
        row << padRight(tcol.str(), kTimeW);
        row << padRight(p.endMem.currentBytes ? formatBytes(p.endMem.currentBytes) : "-", kMemW);
        row << padRight(p.endMem.peakBytes ? formatBytes(p.endMem.peakBytes) : "-", kMemW);
        writeLine(row.str());
    }
    writeLine("  " + std::string(kPhaseW + kTimeW + 2 * kMemW - 2, '-'));

    std::ostringstream totalRow;
    totalRow << "  " << padRight("TOTAL", kPhaseW);
    std::ostringstream tcol;
    tcol << std::fixed << std::setprecision(3) << totalSec << " s";
    totalRow << padRight(tcol.str(), kTimeW);
    totalRow << padRight(finalMem.currentBytes ? formatBytes(finalMem.currentBytes) : "-", kMemW);
    totalRow << padRight(finalMem.peakBytes ? formatBytes(finalMem.peakBytes) : "-", kMemW);
    writeLine(totalRow.str());
}

bool Logger::writeCrashDump(const std::string& crashLogPath, const std::string& whatMessage) {
    std::ofstream crash(crashLogPath, std::ios::out | std::ios::trunc);
    if (!crash.good()) {
        return false;
    }
    crash << "============================================================\n";
    crash << "  CRASH DUMP\n";
    crash << "============================================================\n";
    crash << "what: " << whatMessage << "\n";
    crash << "time: " << timestamp() << "\n";

    // Active phase (the one not yet closed) tells us *where* we were.
    if (!phases_.empty()) {
        const auto& cur = phases_.back();
        crash << "phase: " << cur.name << "\n";
        if (!cur.closed) {
            const double sec = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - cur.start).count();
            crash << "phase_elapsed_s: " << std::fixed << std::setprecision(3) << sec << "\n";
        }
    }
    // Context stack tells us what sub-step (e.g. which frequency).
    if (!context_.empty()) {
        crash << "context:\n";
        for (const auto& [k, v] : context_) {
            crash << "  " << k << " = " << v << "\n";
        }
    }
    const auto mem = queryProcessMemory();
    if (mem.currentBytes > 0) {
        crash << "memory_current: " << formatBytes(mem.currentBytes) << "\n";
        crash << "memory_peak: " << formatBytes(mem.peakBytes) << "\n";
    }
    crash << "\n--- last " << recentLines_.size() << " log lines ---\n";
    for (const auto& l : recentLines_) {
        crash << l << "\n";
    }
    crash.flush();
    return true;
}

void Logger::closeCurrentPhase() {
    if (phases_.empty()) {
        return;
    }
    auto& cur = phases_.back();
    if (cur.closed) {
        return;
    }
    const auto end = std::chrono::steady_clock::now();
    cur.elapsedSec = std::chrono::duration<double>(end - cur.start).count();
    cur.endMem = queryProcessMemory();
    cur.closed = true;
}

void Logger::writeLine(const std::string& line) {
    console_ << line << '\n';
    console_.flush();
    if (file_ && file_->good()) {
        (*file_) << line << '\n';
        file_->flush();
    }
    // Maintain rolling buffer for crash dumps. Keep the last
    // kRecentCapacity lines; drop oldest.
    recentLines_.push_back(line);
    if (recentLines_.size() > kRecentCapacity) {
        recentLines_.pop_front();
    }
}

std::string Logger::timestamp() {
    const auto now = std::chrono::system_clock::now();
    const auto t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#ifdef _WIN32
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream out;
    out << std::put_time(&tm, "%H:%M:%S");
    return out.str();
}

}  // namespace fem
