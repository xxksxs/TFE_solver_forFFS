#include "bpfem/sweep/AlpsSweep.hpp"

#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fem::sweep {

AlpsSweep::AlpsSweep(const ProjectDefinition& project,
                     const FEMAssembler& assembler,
                     const PortModeSolver& portModeSolver,
                     AlpsOptions options)
    : project_(project), assembler_(assembler), portModeSolver_(portModeSolver), options_(options) {}

SweepResult AlpsSweep::run(const std::vector<double>& frequencies, const SweepContext& ctx) {
    const double expansionHz = options_.expansionFrequencyHz > 0.0
        ? options_.expansionFrequencyHz
        : (frequencies.empty()
            ? 0.5 * (project_.sweep.startHz + project_.sweep.endHz)
            : 0.5 * (frequencies.front() + frequencies.back()));

    ctx.log.info("ALPS expansion frequency: " + std::to_string(expansionHz / 1.0e9) + " GHz");
    ctx.log.info("ALPS-Galerkin MVP Krylov order per port: " + std::to_string(options_.krylovOrder));

    linalg::SolverConfig solverCfg;
    solverCfg.maxIterations = ctx.linearMaxIterations;
    solverCfg.tolerance = ctx.linearTolerance;
    const int romDim = buildOffline(expansionHz, ctx.solver, solverCfg);
    ctx.log.info("ALPS ROM dimension: " + std::to_string(romDim));
    ctx.log.info("ALPS retained candidate columns: " + std::to_string(retainedColumns_));
    ctx.log.info("ALPS deflated candidate columns: " + std::to_string(deflatedColumns_));

    SweepResult out;
    out.points.reserve(frequencies.size());
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        const double f = frequencies[i];
        if (i == 0 || (i + 1) % 50 == 0 || i + 1 == frequencies.size()) {
            ctx.log.info("ALPS evaluating point " + std::to_string(i + 1) + "/"
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
        diag.expansionFrequenciesHz = {expansionHz};
        diag.requestedOrder = options_.krylovOrder;
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

int AlpsSweep::buildOffline(double expansionFrequencyHz,
                            linalg::ISparseSolver& solver,
                            const linalg::SolverConfig& solverConfig) {
    expansionFrequencyHz_ = expansionFrequencyHz;
    numProjectPorts_ = static_cast<int>(project_.ports.size());
    if (numProjectPorts_ == 0) {
        throw std::runtime_error("AlpsSweep: no ports defined");
    }

    affine_ = assembler_.buildAffineSystem();
    fullDim_ = static_cast<int>(affine_.K.size());
    numVirtualPorts_ = static_cast<int>(affine_.portCoupling.size());
    if (numVirtualPorts_ == 0) {
        throw std::runtime_error("AlpsSweep: AffineSystem produced zero port entries");
    }
    if (!affine_.lossless) {
        throw std::runtime_error("AlpsSweep: lossy materials (sigma != 0) are not supported in this MVP. Use --sweep direct.");
    }

    std::vector<Complex> rhsTemp;
    SparseMatrix A0 = assembler_.assemble(expansionFrequencyHz_, rhsTemp);
    auto solveAt = [&](const std::vector<Complex>& b) -> std::vector<Complex> {
        SolveResult r = solver.solve(A0, b, solverConfig);
        return r.field;
    };

    const auto portVectors =
        fastsweep::PolynomialPortMomentBuilder::buildPortVectors(affine_, affine_.K.size());

    std::vector<std::vector<Complex>> basis;
    const int targetCols = std::max(1, options_.krylovOrder) * numVirtualPorts_;
    basis.reserve(static_cast<std::size_t>(targetCols));
    retainedColumns_ = 0;
    deflatedColumns_ = 0;

    auto appendIfIndependent = [&](std::vector<Complex>&& w) {
        const double pre = fastsweep::norm2(w);
        if (pre < options_.dropTolerance) {
            ++deflatedColumns_;
            return;
        }
        const double post = fastsweep::modifiedGramSchmidt(w, basis);
        if (post < options_.dropTolerance * std::max(pre, 1.0)) {
            ++deflatedColumns_;
            return;
        }
        for (auto& z : w) {
            z /= post;
        }
        basis.push_back(std::move(w));
        ++retainedColumns_;
    };

    for (int v = 0; v < numVirtualPorts_; ++v) {
        std::vector<Complex> w = solveAt(portVectors[static_cast<std::size_t>(v)]);
        appendIfIndependent(std::move(w));
    }

    std::size_t cursor = 0;
    while (static_cast<int>(basis.size()) < targetCols && cursor < basis.size()) {
        const auto Mv = affine_.M.multiply(basis[cursor]);
        std::vector<Complex> w = solveAt(Mv);
        appendIfIndependent(std::move(w));
        ++cursor;
    }

    romDim_ = static_cast<int>(basis.size());
    if (romDim_ == 0) {
        throw std::runtime_error("AlpsSweep: Krylov subspace collapsed to zero dimension");
    }

    model_.build(project_, affine_, portModeSolver_, std::move(basis), portVectors);
    ready_ = true;
    return romDim_;
}

SParameterPoint AlpsSweep::evaluate(double frequencyHz) const {
    return model_.evaluate(frequencyHz);
}

std::vector<AlpsSweep::Complex> AlpsSweep::reconstructField(double frequencyHz) const {
    return model_.reconstructField(frequencyHz);
}

}  // namespace fem::sweep
