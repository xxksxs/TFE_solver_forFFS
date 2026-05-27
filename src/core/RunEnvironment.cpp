#include "bpfem/core/RunEnvironment.hpp"

#include "bpfem/core/Logger.hpp"

#include <array>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <sstream>
#include <thread>

#if defined(_WIN32)
  #define WIN32_LEAN_AND_MEAN
  #define NOMINMAX
  #include <windows.h>
  #include <intrin.h>
#elif defined(__linux__)
  #include <sys/utsname.h>
  #include <unistd.h>
  #include <fstream>
#endif

#ifdef BPFEM_USE_MKL
  #include <mkl_service.h>
#endif

namespace fem {

namespace {

// CPU brand string via CPUID leaves 0x80000002 - 0x80000004 on x86; falls
// back to "n/a" elsewhere.
std::string queryCpuBrand() {
#if defined(_WIN32) || defined(__GNUC__)
  #if defined(_MSC_VER)
    int cpuInfo[4] = {-1};
    __cpuid(cpuInfo, 0x80000000);
    const unsigned maxExt = static_cast<unsigned>(cpuInfo[0]);
    if (maxExt < 0x80000004U) {
        return "n/a";
    }
    char brand[0x40] = {};
    for (unsigned i = 0; i < 3; ++i) {
        __cpuid(cpuInfo, static_cast<int>(0x80000002U + i));
        std::memcpy(brand + i * 16, cpuInfo, sizeof(cpuInfo));
    }
    // The brand string is space-padded; trim leading whitespace.
    std::string s(brand);
    const auto first = s.find_first_not_of(" \t");
    return (first == std::string::npos) ? std::string("n/a") : s.substr(first);
  #else
    return "n/a";  // GCC/Clang on x86 could use __get_cpuid; out of MVP scope
  #endif
#else
    return "n/a";
#endif
}

std::string queryHostName() {
#if defined(_WIN32)
    char buf[MAX_COMPUTERNAME_LENGTH + 1] = {};
    DWORD sz = sizeof(buf);
    if (GetComputerNameA(buf, &sz)) {
        return std::string(buf);
    }
#elif defined(__linux__)
    char buf[256] = {};
    if (gethostname(buf, sizeof(buf) - 1) == 0) {
        return std::string(buf);
    }
#endif
    return "n/a";
}

std::string queryOsDescription() {
#if defined(_WIN32)
    OSVERSIONINFOEXA info{};
    info.dwOSVersionInfoSize = sizeof(info);
    // GetVersionExA is deprecated and lies about >= Win 8.1, but for log-
    // header purposes the major.minor + build is enough; we explicitly
    // ignore the deprecation warning here.
    #pragma warning(push)
    #pragma warning(disable: 4996)
    if (GetVersionExA(reinterpret_cast<OSVERSIONINFOA*>(&info))) {
        std::ostringstream out;
        out << "Windows " << info.dwMajorVersion << "." << info.dwMinorVersion
            << " build " << info.dwBuildNumber;
        return out.str();
    }
    #pragma warning(pop)
    return "Windows (unknown version)";
#elif defined(__linux__)
    utsname u{};
    if (uname(&u) == 0) {
        std::ostringstream out;
        out << u.sysname << " " << u.release;
        return out.str();
    }
    return "Linux (unknown)";
#elif defined(__APPLE__)
    return "macOS";
#else
    return "n/a";
#endif
}

std::string queryMemoryDescription() {
#if defined(_WIN32)
    MEMORYSTATUSEX info{};
    info.dwLength = sizeof(info);
    if (GlobalMemoryStatusEx(&info)) {
        std::ostringstream out;
        const double totalGb = info.ullTotalPhys / (1024.0 * 1024.0 * 1024.0);
        const double availGb = info.ullAvailPhys / (1024.0 * 1024.0 * 1024.0);
        out.setf(std::ios::fixed);
        out.precision(2);
        out << totalGb << " GB total, " << availGb << " GB available";
        return out.str();
    }
#elif defined(__linux__)
    std::ifstream meminfo("/proc/meminfo");
    long long memTotalKb = 0, memAvailKb = 0;
    std::string line;
    while (std::getline(meminfo, line)) {
        if (line.rfind("MemTotal:", 0) == 0) {
            std::sscanf(line.c_str(), "MemTotal: %lld kB", &memTotalKb);
        } else if (line.rfind("MemAvailable:", 0) == 0) {
            std::sscanf(line.c_str(), "MemAvailable: %lld kB", &memAvailKb);
        }
    }
    if (memTotalKb > 0) {
        std::ostringstream out;
        out.setf(std::ios::fixed);
        out.precision(2);
        out << (memTotalKb / (1024.0 * 1024.0)) << " GB total, "
            << (memAvailKb / (1024.0 * 1024.0)) << " GB available";
        return out.str();
    }
#endif
    return "n/a";
}

std::string queryRuntimeDescription() {
    std::ostringstream out;
#ifdef BPFEM_USE_MKL
    char mklVersion[256] = {};
    mkl_get_version_string(mklVersion, sizeof(mklVersion));
    // The MKL string contains trailing spaces; trim.
    std::string mkl(mklVersion);
    while (!mkl.empty() && mkl.back() == ' ') mkl.pop_back();
    out << "MKL: " << mkl;
    out << ", MKL threads=" << mkl_get_max_threads();
#else
    out << "MKL: disabled";
#endif
    out << ", std::thread::hardware_concurrency=" << std::thread::hardware_concurrency();
    return out.str();
}

std::string queryBuildString() {
    std::ostringstream out;
#if defined(_MSC_VER)
    out << "MSVC " << (_MSC_VER / 100) << "." << (_MSC_VER % 100);
#elif defined(__clang__)
    out << "Clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#elif defined(__GNUC__)
    out << "GCC " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#else
    out << "unknown compiler";
#endif

#if defined(NDEBUG)
    out << " / Release";
#else
    out << " / Debug";
#endif

#ifdef BPFEM_USE_MKL
    out << " / BPFEM_USE_MKL=ON";
#else
    out << " / BPFEM_USE_MKL=OFF";
#endif

#ifdef BPFEM_OUTPUT_PARAVIEW
    out << " / BPFEM_OUTPUT_PARAVIEW=ON";
#else
    out << " / BPFEM_OUTPUT_PARAVIEW=OFF";
#endif

    return out.str();
}

std::string joinCommandLine(int argc, char** argv) {
    std::ostringstream out;
    for (int i = 0; i < argc; ++i) {
        if (i > 0) out << ' ';
        const std::string a = argv[i];
        // Quote any argv that contains a space so the log line is
        // copy-paste reproducible.
        if (a.find(' ') != std::string::npos) {
            out << '"' << a << '"';
        } else {
            out << a;
        }
    }
    return out.str();
}

}  // namespace

RunEnvironment captureRunEnvironment(int argc, char** argv) {
    RunEnvironment env;
    env.solverVersion = "bp_fem_solver 2026R1";
    env.buildString = queryBuildString();
    env.runtimeString = queryRuntimeDescription();
    env.hostName = queryHostName();
    env.osDescription = queryOsDescription();
    {
        const std::string brand = queryCpuBrand();
        std::ostringstream cpuOut;
        cpuOut << brand << ", " << std::thread::hardware_concurrency() << " logical cores";
        env.cpuDescription = cpuOut.str();
    }
    env.memoryDescription = queryMemoryDescription();
    env.commandLine = joinCommandLine(argc, argv);
    try {
        env.workingDirectory = std::filesystem::current_path().string();
    } catch (...) {
        env.workingDirectory = "n/a";
    }
    return env;
}

void writeEnvHeader(Logger& log, const RunEnvironment& env) {
    log.env("solver", env.solverVersion);
    log.env("build", env.buildString);
    log.env("runtime", env.runtimeString);
    log.env("host", env.hostName + ", " + env.osDescription);
    log.env("cpu", env.cpuDescription);
    log.env("memory", env.memoryDescription);
    log.env("cmdline", env.commandLine);
    log.env("cwd", env.workingDirectory);
}

}  // namespace fem
