#include "bpfem/linalg/MumpsBackend.hpp"
#include "bpfem/linalg/SolveResidual.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef BPFEM_USE_MUMPS
#include <zmumps_c.h>
#endif

namespace fem::linalg {
namespace {

#ifdef BPFEM_USE_MUMPS
using Clock = std::chrono::steady_clock;

double elapsedSince(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

MUMPS_INT toMumpsInt(std::size_t value) {
    if (value > static_cast<std::size_t>(std::numeric_limits<MUMPS_INT>::max())) {
        throw std::overflow_error("MUMPS dimension or index is too large");
    }
    return static_cast<MUMPS_INT>(value);
}
#endif

}  // namespace

struct MumpsBackend::State {
#ifdef BPFEM_USE_MUMPS
    ZMUMPS_STRUC_C handle{};
    std::vector<MUMPS_INT> rows;
    std::vector<MUMPS_INT> columns;
    std::vector<mumps_double_complex> values;

    void run(MUMPS_INT job) {
        handle.job = job;
        zmumps_c(&handle);
        if (handle.infog[0] < 0) {
            throw std::runtime_error("ZMUMPS job " + std::to_string(job)
                + " failed, INFOG(1)=" + std::to_string(handle.infog[0])
                + ", INFOG(2)=" + std::to_string(handle.infog[1]));
        }
    }

    ~State() {
        handle.job = -2;
        zmumps_c(&handle);
    }
#endif
};

MumpsBackend::MumpsBackend(MumpsMode mode, double blrTolerance, std::filesystem::path oocDirectory)
    : mode_(mode), blrTolerance_(blrTolerance), oocDirectory_(std::move(oocDirectory)) {
#ifndef BPFEM_USE_MUMPS
    throw std::runtime_error("MUMPS requires a build with -DBPFEM_USE_MUMPS=ON");
#endif
}

MumpsBackend::~MumpsBackend() = default;

void MumpsBackend::clearFactorization() {
    state_.reset();
}

void MumpsBackend::factorize(const SparseMatrix& matrix, const SolverConfig&) {
#ifdef BPFEM_USE_MUMPS
    clearFactorization();
    if (matrix.size() == 0) {
        throw std::invalid_argument("MUMPS requires a nonempty matrix");
    }
    auto state = std::make_unique<State>();
    auto& mumps = state->handle;
    mumps.par = 1;
    mumps.sym = 2;
    mumps.comm_fortran = -987654;  // MUMPS sentinel for its default communicator.
    state->run(-1);
    std::fill_n(mumps.icntl, 4, 0);
    // Keep tree-based OpenMP scheduling off. The sequential-MPI build can
    // still use its other OpenMP kernels without reordering the tree for
    // a particular thread count.
    mumps.icntl[47] = 0;
    mumps.icntl[21] = mode_ == MumpsMode::OutOfCore ? 1 : 0;
    if (mode_ == MumpsMode::OutOfCore) {
        const auto directory = oocDirectory_.empty() ? std::filesystem::current_path() : oocDirectory_;
        const auto path = directory.string();
        if (!std::filesystem::is_directory(directory) || path.size() >= sizeof(mumps.ooc_tmpdir)) {
            throw std::invalid_argument("Invalid MUMPS out-of-core directory");
        }
        std::copy(path.c_str(), path.c_str() + path.size() + 1, mumps.ooc_tmpdir);
        const auto prefix = "bpfem-" + std::to_string(Clock::now().time_since_epoch().count());
        std::copy(prefix.c_str(), prefix.c_str() + prefix.size() + 1, mumps.ooc_prefix);
    }
    mumps.icntl[34] = mode_ == MumpsMode::BlockLowRank ? 2 : 0;
    if (mode_ == MumpsMode::BlockLowRank) {
        mumps.cntl[6] = blrTolerance_;
    }

    const auto& offsets = matrix.rowOffsets();
    const auto& columns = matrix.columnIndices();
    const auto& values = matrix.values();
    mumps.n = toMumpsInt(matrix.size());
    mumps.nnz = static_cast<MUMPS_INT8>(values.size());
    state->rows.reserve(values.size());
    state->columns.reserve(values.size());
    state->values.reserve(values.size());
    for (std::size_t row = 0; row < matrix.size(); ++row) {
        for (int i = offsets[row]; i < offsets[row + 1]; ++i) {
            state->rows.push_back(toMumpsInt(static_cast<std::size_t>(columns[i]) + 1));
            state->columns.push_back(toMumpsInt(row + 1));
            state->values.push_back({values[i].real(), values[i].imag()});
        }
    }
    mumps.irn = state->rows.data();
    mumps.jcn = state->columns.data();
    mumps.a = state->values.data();
    auto started = Clock::now();
    state->run(1);
    statistics_.symbolicAnalysisSec += elapsedSince(started);
    ++statistics_.symbolicAnalysisCount;
    started = Clock::now();
    state->run(2);
    statistics_.numericFactorizationSec += elapsedSince(started);
    ++statistics_.numericFactorizationCount;
    state_ = std::move(state);
#else
    (void)matrix;
    throw std::logic_error("MUMPS support unavailable");
#endif
}

SolveResult MumpsBackend::solve(
    const SparseMatrix& matrix, const std::vector<std::complex<double>>& rhs,
    const SolverConfig& config) {
#ifdef BPFEM_USE_MUMPS
    factorize(matrix, config);
    auto result = solveFactorized(rhs, config);
    checkResidual(matrix, rhs, config.tolerance, result);
    return result;
#else
    (void)matrix;
    (void)rhs;
    (void)config;
    throw std::logic_error("MUMPS support unavailable");
#endif
}

SolveResult MumpsBackend::solveFactorized(
    const std::vector<std::complex<double>>& rhs, const SolverConfig& config) {
    auto results = solveFactorizedBatch({rhs}, config);
    return std::move(results.front());
}

std::vector<SolveResult> MumpsBackend::solveFactorizedBatch(
    const std::vector<std::vector<std::complex<double>>>& rhsList, const SolverConfig&) {
#ifdef BPFEM_USE_MUMPS
    if (!state_) {
        throw std::logic_error("MUMPS solve requires a successful factorization");
    }
    if (rhsList.empty()) {
        return {};
    }
    auto& mumps = state_->handle;
    const auto n = static_cast<std::size_t>(mumps.n);
    const auto nrhs = toMumpsInt(rhsList.size());
    std::vector<mumps_double_complex> values;
    if (rhsList.size() > values.max_size() / n) {
        throw std::overflow_error("MUMPS RHS batch is too large");
    }
    values.reserve(n * rhsList.size());
    for (const auto& rhs : rhsList) {
        if (rhs.size() != n) {
            throw std::invalid_argument("MUMPS RHS size does not match the matrix");
        }
        for (const auto& value : rhs) {
            values.push_back({value.real(), value.imag()});
        }
    }
    mumps.nrhs = nrhs;
    mumps.lrhs = mumps.n;
    mumps.rhs = values.data();
    const auto started = Clock::now();
    try {
        state_->run(3);
    } catch (...) {
        mumps.rhs = nullptr;
        clearFactorization();
        throw;
    }
    mumps.rhs = nullptr;
    statistics_.factorizedRhsSolveSec += elapsedSince(started);
    ++statistics_.factorizedSolveCallCount;
    statistics_.factorizedRhsSolveCount += rhsList.size();
    statistics_.batchRhsMax = std::max<std::uint64_t>(statistics_.batchRhsMax, rhsList.size());
    std::vector<SolveResult> results(rhsList.size());
    for (std::size_t k = 0; k < results.size(); ++k) {
        auto& result = results[k];
        result.iterations = 1;
        result.field.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& value = values[k * n + i];
            result.field.emplace_back(value.r, value.i);
        }
    }
    return results;
#else
    (void)rhsList;
    throw std::logic_error("MUMPS support unavailable");
#endif
}

}  // namespace fem::linalg
