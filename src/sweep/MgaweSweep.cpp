#include "bpfem/sweep/MgaweSweep.hpp"

#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"
#include "bpfem/fastsweep/PortModeUtilities.hpp"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fem::sweep {

MgaweSweep::MgaweSweep(const ProjectDefinition& project,
                       const FEMAssembler& assembler,
                       const PortModeSolver& portModeSolver,
                       MgaweOptions options)
    : project_(project),
      assembler_(assembler),
      portModeSolver_(portModeSolver),
      options_(options) {
    options_.expansionPointCount = std::max(1, options_.expansionPointCount);
    options_.localOrder = std::max(1, options_.localOrder);
    if (options_.dropTolerance <= 0.0) {
        options_.dropTolerance = 1.0e-10;
    }
}

SweepResult MgaweSweep::run(const std::vector<double>& frequencies, const SweepContext& ctx) {
    if (frequencies.empty()) {
        return {};
    }

    const double left = frequencies.front();
    const double right = frequencies.back();
    const int pointCount = std::max(1, options_.expansionPointCount);
    std::vector<double> expansionFrequencies;
    expansionFrequencies.reserve(static_cast<std::size_t>(pointCount));
    if (pointCount == 1 || left == right) {
        expansionFrequencies.push_back(0.5 * (left + right));
    } else {
        for (int i = 0; i < pointCount; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(pointCount - 1);
            expansionFrequencies.push_back(left + t * (right - left));
        }
    }

    std::ostringstream points;
    points << std::fixed << std::setprecision(6);
    for (std::size_t i = 0; i < expansionFrequencies.size(); ++i) {
        if (i != 0) {
            points << ", ";
        }
        points << (expansionFrequencies[i] / 1.0e9);
    }
    ctx.log.info("MGAWE expansion frequencies: [" + points.str() + "] GHz");
    ctx.log.info("MGAWE local AWE moment count per point: " + std::to_string(options_.localOrder));
    std::ostringstream dropTol;
    dropTol << std::scientific << std::setprecision(3) << options_.dropTolerance;
    ctx.log.info("MGAWE global MGS drop tolerance: " + dropTol.str());

    linalg::SolverConfig solverCfg;
    solverCfg.maxIterations = ctx.linearMaxIterations;
    solverCfg.tolerance = ctx.linearTolerance;
    const int dim = buildOffline(expansionFrequencies, ctx.solver, solverCfg);
    ctx.log.info("MGAWE ROM dimension: " + std::to_string(dim));
    ctx.log.info("MGAWE retained candidate columns: " + std::to_string(retainedColumns_));
    ctx.log.info("MGAWE deflated candidate columns: " + std::to_string(deflatedColumns_));

    SweepResult out;
    out.points.reserve(frequencies.size());
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        const double f = frequencies[i];
        if (i == 0 || (i + 1) % 50 == 0 || i + 1 == frequencies.size()) {
            ctx.log.info("MGAWE evaluating point " + std::to_string(i + 1) + "/"
                         + std::to_string(frequencies.size())
                         + " @ " + std::to_string(f / 1.0e9) + " GHz");
        }
        out.points.push_back(evaluate(f));
        out.lastFrequencyHz = f;
    }
    if (!frequencies.empty()) {
        out.lastEdgeDofs = reconstructField(frequencies.back());
    }
    if (!ctx.outputDirectory.empty()) {
        fastsweep::FastSweepDiagnostics diag;
        diag.algorithm = name();
        diag.expansionFrequenciesHz = expansionFrequencies;
        diag.requestedOrder = options_.localOrder;
        diag.romDimension = romDim_;
        diag.retainedColumns = retainedColumns_;
        diag.deflatedColumns = deflatedColumns_;
        diag.basisOrthogonalityError = model_.basisOrthogonalityError();
        diag.reducedSolveSucceeded = model_.lastSolveSucceeded();
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

int MgaweSweep::buildOffline(const std::vector<double>& expansionFrequencies,
                             linalg::ISparseSolver& solver,
                             const linalg::SolverConfig& solverConfig) {
    if (project_.ports.size() < 2) {
        throw std::runtime_error("MgaweSweep: at least two project ports are required");
    }
    if (expansionFrequencies.empty()) {
        throw std::runtime_error("MgaweSweep: at least one expansion point is required");
    }

    numProjectPorts_ = static_cast<int>(project_.ports.size());
    affine_ = assembler_.buildAffineSystem();
    fullDim_ = static_cast<int>(affine_.K.size());
    numVirtualPorts_ = static_cast<int>(affine_.portCoupling.size());
    if (numVirtualPorts_ == 0) {
        throw std::runtime_error("MgaweSweep: AffineSystem produced zero port entries");
    }
    if (!affine_.lossless) {
        throw std::runtime_error("MgaweSweep: lossy materials (sigma != 0) are not supported. Use --sweep direct.");
    }

    const auto portVectors =
        fastsweep::PolynomialPortMomentBuilder::buildPortVectors(affine_, affine_.K.size());

    basis_.clear();
    retainedColumns_ = 0;
    deflatedColumns_ = 0;
    const int candidateCount =
        static_cast<int>(expansionFrequencies.size()) * std::max(1, options_.localOrder);
    basis_.reserve(static_cast<std::size_t>(candidateCount));

    for (double expansionHz : expansionFrequencies) {
        const auto localMoments = buildLocalMoments(expansionHz, solver, solverConfig, portVectors);
        for (const auto& moment : localMoments) {
            appendIfIndependent(std::vector<Complex>(moment.begin(), moment.end()));
        }
    }

    romDim_ = static_cast<int>(basis_.size());
    if (romDim_ == 0) {
        throw std::runtime_error("MgaweSweep: global basis collapsed to zero dimension");
    }

    projectReducedModel(portVectors);
    ready_ = true;
    return romDim_;
}

std::vector<std::vector<MgaweSweep::Complex>> MgaweSweep::buildLocalMoments(
    double expansionFrequencyHz,
    linalg::ISparseSolver& solver,
    const linalg::SolverConfig& solverConfig,
    const std::vector<std::vector<Complex>>& portVectors) const {
    return fastsweep::PolynomialPortMomentBuilder::generateLosslessMoments(
        project_, assembler_, portModeSolver_, affine_, expansionFrequencyHz,
        std::max(1, options_.localOrder), solver, solverConfig, portVectors);
}

bool MgaweSweep::appendIfIndependent(std::vector<Complex>&& w) {
    const double pre = fastsweep::norm2(w);
    if (pre < options_.dropTolerance) {
        ++deflatedColumns_;
        return false;
    }
    const double post = fastsweep::modifiedGramSchmidt(w, basis_);
    if (post < options_.dropTolerance * std::max(pre, 1.0)) {
        ++deflatedColumns_;
        return false;
    }
    for (auto& z : w) {
        z /= post;
    }
    basis_.push_back(std::move(w));
    ++retainedColumns_;
    return true;
}

void MgaweSweep::projectReducedModel(const std::vector<std::vector<Complex>>& portVectors) {
    model_.build(project_, affine_, portModeSolver_, std::move(basis_), portVectors);
}

std::vector<MgaweSweep::Complex> MgaweSweep::solveReduced(double frequencyHz) const {
    return model_.solveReduced(frequencyHz);
}

SParameterPoint MgaweSweep::evaluate(double frequencyHz) const {
    return model_.evaluate(frequencyHz);
}

std::vector<MgaweSweep::Complex> MgaweSweep::reconstructField(double frequencyHz) const {
    return model_.reconstructField(frequencyHz);
}

const PortMode& MgaweSweep::virtualPortMode(int virtualPortIndex) const {
    return fastsweep::virtualPortMode(portModeSolver_, affine_, virtualPortIndex);
}

}  // namespace fem::sweep
