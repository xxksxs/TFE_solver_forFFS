#include "bpfem/sweep/AweSweep.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"
#include "bpfem/fastsweep/PortModeUtilities.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem::sweep {

namespace {

using Complex = std::complex<double>;

}  // namespace

AweSweep::AweSweep(const ProjectDefinition& project,
                   const FEMAssembler& assembler,
                   const PortModeSolver& portModeSolver,
                   AweOptions options)
    : project_(project), assembler_(assembler), portModeSolver_(portModeSolver), options_(options) {
    if (options_.order < 1) {
        options_.order = 1;
    }
}

SweepResult AweSweep::run(const std::vector<double>& frequencies, const SweepContext& ctx) {
    const double expansionHz = options_.expansionFrequencyHz > 0.0
        ? options_.expansionFrequencyHz
        : (frequencies.empty()
            ? 0.5 * (project_.sweep.startHz + project_.sweep.endHz)
            : 0.5 * (frequencies.front() + frequencies.back()));

    ctx.log.info("AWE expansion frequency: " + std::to_string(expansionHz / 1.0e9) + " GHz");
    ctx.log.info("AWE Padé denominator order: " + std::to_string(options_.order));
    ctx.log.info("AWE moment variable: normalized lambda shift (lambda-lambda0)/lambda0");

    linalg::SolverConfig solverCfg;
    solverCfg.maxIterations = ctx.linearMaxIterations;
    solverCfg.tolerance = ctx.linearTolerance;
    buildOffline(expansionHz, ctx.solver, solverCfg);

    ctx.log.info("AWE input Padé order: ["
                 + std::to_string(inputProjectionPade_.numeratorOrder()) + "/"
                 + std::to_string(inputProjectionPade_.denominatorOrder()) + "]");
    ctx.log.info("AWE output Padé order: ["
                 + std::to_string(outputProjectionPade_.numeratorOrder()) + "/"
                 + std::to_string(outputProjectionPade_.denominatorOrder()) + "]");

    SweepResult out;
    out.points.reserve(frequencies.size());
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        const double f = frequencies[i];
        if (i == 0 || (i + 1) % 50 == 0 || i + 1 == frequencies.size()) {
            ctx.log.info("AWE evaluating point " + std::to_string(i + 1) + "/"
                         + std::to_string(frequencies.size())
                         + " @ " + std::to_string(f / 1.0e9) + " GHz");
        }
        out.points.push_back(evaluate(f));
        out.lastFrequencyHz = f;
    }
    if (!ctx.outputDirectory.empty()) {
        fastsweep::FastSweepDiagnostics diag;
        diag.algorithm = name();
        diag.expansionFrequenciesHz = {expansionHz};
        diag.requestedOrder = options_.order;
        diag.padeInputPivotRatio = inputProjectionPade_.pivotRatio();
        diag.padeOutputPivotRatio = outputProjectionPade_.pivotRatio();
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

int AweSweep::buildOffline(double expansionFrequencyHz,
                           linalg::ISparseSolver& solver,
                           const linalg::SolverConfig& solverConfig) {
    if (project_.ports.size() < 2) {
        throw std::runtime_error("AweSweep: at least two project ports are required");
    }

    expansionFrequencyHz_ = expansionFrequencyHz;
    const double kExpansion = 2.0 * pi * expansionFrequencyHz_ / c0;
    lambda0_ = kExpansion * kExpansion;
    lambdaScale_ = std::max(lambda0_, 1.0);

    affine_ = assembler_.buildAffineSystem();
    if (!affine_.lossless) {
        throw std::runtime_error("AweSweep: lossy materials (sigma != 0) are not supported in this MVP. Use --sweep direct.");
    }

    const int q = std::max(1, options_.order);
    const int momentCount = 2 * q;
    const auto portVectors =
        fastsweep::PolynomialPortMomentBuilder::buildPortVectors(affine_, affine_.K.size());
    const auto fieldMoments =
        fastsweep::PolynomialPortMomentBuilder::generateLosslessMoments(
            project_, assembler_, portModeSolver_, affine_, expansionFrequencyHz_,
            momentCount, solver, solverConfig, portVectors);

    std::vector<Complex> inputProjectionMoments;
    std::vector<Complex> outputProjectionMoments;
    inputProjectionMoments.reserve(fieldMoments.size());
    outputProjectionMoments.reserve(fieldMoments.size());
    const int inputFaceId = project_.ports[0].faceId;
    const int outputFaceId = project_.ports[1].faceId;
    for (const auto& xj : fieldMoments) {
        inputProjectionMoments.push_back(portProjection(xj, inputFaceId));
        outputProjectionMoments.push_back(portProjection(xj, outputFaceId));
    }

    inputProjectionPade_ = fastsweep::PadeApproximant::buildBestEffort(
        inputProjectionMoments, q - 1, q);
    outputProjectionPade_ = fastsweep::PadeApproximant::buildBestEffort(
        outputProjectionMoments, q - 1, q);

    ready_ = true;
    return momentCount;
}

SParameterPoint AweSweep::evaluate(double frequencyHz) const {
    if (!ready_) {
        throw std::runtime_error("AweSweep::evaluate called before buildOffline");
    }

    SParameterPoint sp;
    sp.frequencyHz = frequencyHz;
    if (project_.ports.size() < 2) {
        return sp;
    }

    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double lambda = k0 * k0;
    const Complex t((lambda - lambda0_) / lambdaScale_, 0.0);

    const PortMode& inputMode = portModeSolver_.solve(project_.ports[0].faceId);
    const PortMode& outputMode = portModeSolver_.solve(project_.ports[1].faceId);
    const double sIn = powerNormalizationFactor(inputMode, frequencyHz);
    const double sOut = powerNormalizationFactor(outputMode, frequencyHz);
    if (sIn <= 0.0) {
        return sp;
    }

    const Complex rawInput = inputProjectionPade_.evaluate(t);
    const Complex rawOutput = outputProjectionPade_.evaluate(t);
    const Complex bIn = rawInput / sIn;
    const Complex bOut = (sOut > 0.0) ? rawOutput / sOut : Complex(0.0, 0.0);

    const auto& port0 = project_.ports[0];
    const Complex incident = port0.excited
        ? std::polar(std::sqrt(std::max(port0.magnitudeW, 0.0)), port0.phaseDeg * pi / 180.0)
        : Complex(1.0, 0.0);
    if (incident == Complex(0.0, 0.0)) {
        return sp;
    }

    sp.s11 = (bIn - incident) / incident;
    sp.s21 = bOut / incident;
    return sp;
}

AweSweep::Complex AweSweep::portProjection(const std::vector<Complex>& edgeDofs, int faceId) const {
    const auto& mode = portModeSolver_.solve(faceId).couplingWeights;
    Complex sum(0.0, 0.0);
    for (const auto& [edgeIndex, value] : mode) {
        if (edgeIndex >= 0 && static_cast<std::size_t>(edgeIndex) < edgeDofs.size()) {
            sum += edgeDofs[static_cast<std::size_t>(edgeIndex)] * value;
        }
    }
    return sum;
}

const PortMode& AweSweep::virtualPortMode(int virtualPortIndex) const {
    return fastsweep::virtualPortMode(portModeSolver_, affine_, virtualPortIndex);
}

}  // namespace fem::sweep
