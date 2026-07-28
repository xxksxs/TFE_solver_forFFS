#include "bpfem/sweep/GaweSweep.hpp"

#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace fem::sweep {

namespace {

GaweOptions normalize(GaweOptions options) {
    options.order = std::max(1, options.order);
    if (options.dropTolerance <= 0.0) {
        options.dropTolerance = 1.0e-10;
    }
    return options;
}

MgaweOptions toSinglePointMgawe(const GaweOptions& options) {
    MgaweOptions mgawe;
    mgawe.expansionPointCount = 1;
    mgawe.localOrder = options.order;
    mgawe.dropTolerance = options.dropTolerance;
    return mgawe;
}

}  // namespace

GaweSweep::GaweSweep(const ProjectDefinition& project,
                     const FEMAssembler& assembler,
                     const PortModeSolver& portModeSolver,
                     GaweOptions options)
    : options_(normalize(options)),
      delegate_(project, assembler, portModeSolver, toSinglePointMgawe(options_)) {}

SweepResult GaweSweep::run(const std::vector<double>& frequencies, const SweepContext& ctx) {
    const double expansionHz = options_.expansionFrequencyHz > 0.0
        ? options_.expansionFrequencyHz
        : (frequencies.empty()
            ? 0.5 * (ctx.project.sweep.startHz + ctx.project.sweep.endHz)
            : 0.5 * (frequencies.front() + frequencies.back()));

    ctx.log.info("GAWE expansion frequency: " + std::to_string(expansionHz / 1.0e9) + " GHz");
    ctx.log.info("GAWE local AWE/Galerkin order: " + std::to_string(options_.order));
    std::ostringstream dropTol;
    dropTol << std::scientific << std::setprecision(3) << options_.dropTolerance;
    ctx.log.info("GAWE MGS drop tolerance: " + dropTol.str());

    linalg::SolverConfig solverCfg;
    solverCfg.maxIterations = ctx.linearMaxIterations;
    solverCfg.tolerance = ctx.linearTolerance;
    const auto offlineStarted = std::chrono::steady_clock::now();
    const int dim = buildOffline(expansionHz, ctx.solver, solverCfg);
    const double offlineBuildSec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - offlineStarted).count();
    ctx.log.info("GAWE ROM dimension: " + std::to_string(dim));
    ctx.log.info("GAWE retained candidate columns: " + std::to_string(retainedColumns()));
    ctx.log.info("GAWE deflated candidate columns: " + std::to_string(deflatedColumns()));

    SweepResult out;
    out.offlineBuildSec = offlineBuildSec;
    out.orthogonalizationSec = delegate_.lastOrthogonalizationSec();
    out.romProjectionSec = delegate_.lastRomProjectionSec();
    out.points.reserve(frequencies.size());
    const auto onlineStarted = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        const double f = frequencies[i];
        if (i == 0 || (i + 1) % 50 == 0 || i + 1 == frequencies.size()) {
            ctx.log.info("GAWE evaluating point " + std::to_string(i + 1) + "/"
                         + std::to_string(frequencies.size())
                         + " @ " + std::to_string(f / 1.0e9) + " GHz");
        }
        out.points.push_back(evaluate(f));
        out.lastFrequencyHz = f;
    }
    out.onlineSweepSec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - onlineStarted).count();
    if (!frequencies.empty()) {
        out.lastEdgeDofs = reconstructField(frequencies.back());
    }
    if (!ctx.outputDirectory.empty()) {
        fastsweep::FastSweepDiagnostics diag;
        diag.algorithm = name();
        diag.expansionFrequenciesHz = {expansionHz};
        diag.requestedOrder = options_.order;
        diag.romDimension = dimension();
        diag.retainedColumns = retainedColumns();
        diag.deflatedColumns = deflatedColumns();
        diag.basisOrthogonalityError = delegate_.basisOrthogonalityError();
        diag.reducedSolveSucceeded = delegate_.reducedSolveSucceeded();
        diag.maxPassivityError = fastsweep::maxPassivityError(out.points);
        const auto path = ctx.outputDirectory / "diagnostics.json";
        if (fastsweep::writeDiagnosticsJson(path, diag)) {
            ctx.log.info("Wrote " + path.string());
        } else {
            ctx.log.warn("Failed to write " + path.string());
        }
    }
    (void)ctx.onFieldSolved;
    return out;
}

int GaweSweep::buildOffline(double expansionFrequencyHz,
                            linalg::ISparseSolver& solver,
                            const linalg::SolverConfig& solverConfig) {
    return delegate_.buildOffline({expansionFrequencyHz}, solver, solverConfig);
}

SParameterPoint GaweSweep::evaluate(double frequencyHz) const {
    return delegate_.evaluate(frequencyHz);
}

std::vector<GaweSweep::Complex> GaweSweep::reconstructField(double frequencyHz) const {
    return delegate_.reconstructField(frequencyHz);
}

}  // namespace fem::sweep
