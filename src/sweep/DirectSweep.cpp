#include "bpfem/sweep/DirectSweep.hpp"

#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"

#include <chrono>
#include <complex>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace fem::sweep {

namespace {

// Format wall-clock seconds as a fixed-precision string for log lines.
std::string formatSec(double s) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << s << " s";
    return oss.str();
}

}  // namespace

SweepResult DirectSweep::run(const std::vector<double>& frequencies,
                             const SweepContext& ctx) {
    SweepResult out;
    out.points.reserve(frequencies.size());

    linalg::SolverConfig solverCfg;
    solverCfg.maxIterations = ctx.linearMaxIterations;
    solverCfg.tolerance = ctx.linearTolerance;

    using clk = std::chrono::steady_clock;
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        const double f = frequencies[i];

        // Push a context tag so any crash-dump in the subsequent assemble /
        // solve labels which point we were on.
        std::ostringstream tag;
        tag << (i + 1) << "/" << frequencies.size() << " @ " << (f / 1.0e9) << " GHz";
        ctx.log.pushContext("freq", tag.str());

        ctx.log.info("Solving point " + tag.str());

        // 1. Assemble: build A(f) and rhs in CSR.
        const auto t0 = clk::now();
        std::vector<std::complex<double>> rhs;
        SparseMatrix matrix = ctx.assembler.assemble(f, rhs);
        const auto t1 = clk::now();
        const double assembleSec = std::chrono::duration<double>(t1 - t0).count();
        ctx.log.info("  assemble: " + formatSec(assembleSec) + " ("
                     + std::to_string(matrix.values().size()) + " nnz upper-tri)");

        // 2. Solve: PARDISO numeric factor + back-substitute (or iterative).
        ctx.log.info("  solving (this may take seconds for direct factorization)...");
        const auto t2 = clk::now();
        SolveResult result = ctx.solver.solve(matrix, rhs, solverCfg);
        const auto t3 = clk::now();
        const double solveSec = std::chrono::duration<double>(t3 - t2).count();
        ctx.log.info("  solve:    " + formatSec(solveSec)
                     + ", iterations=" + std::to_string(result.iterations)
                     + ", residual=" + std::to_string(result.residual));

        // 3. Extract S parameters from the solved field.
        out.points.push_back(ctx.extractor.extract(f, result.field));
        if (ctx.onFrequencyStat) {
            ctx.onFrequencyStat(f, assembleSec, solveSec,
                                result.iterations, result.residual);
        }
        if (ctx.onFieldSolved) {
            ctx.onFieldSolved(f, result.field);
        }
        out.lastEdgeDofs = std::move(result.field);
        out.lastFrequencyHz = f;

        ctx.log.popContext();
    }

    if (!ctx.outputDirectory.empty()) {
        fastsweep::FastSweepDiagnostics diag;
        diag.algorithm = name();
        diag.romDimension = static_cast<int>(out.lastEdgeDofs.size());
        diag.reducedSolveSucceeded = true;
        diag.maxPassivityError = fastsweep::maxPassivityError(out.points);
        const auto path = ctx.outputDirectory / "diagnostics.json";
        if (fastsweep::writeDiagnosticsJson(path, diag)) {
            ctx.log.info("Wrote " + path.string());
        } else {
            ctx.log.warn("Failed to write " + path.string());
        }
    }
    return out;
}

}  // namespace fem::sweep
