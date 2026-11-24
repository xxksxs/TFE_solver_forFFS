#include "bpfem/sweep/AweSweep.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"
#include "bpfem/fastsweep/PortModeUtilities.hpp"
#include "bpfem/post/ResultExtractor.hpp"

#include <algorithm>
#include <chrono>
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
    ctx.log.info("AWE moment variable: normalized wavenumber shift (k-k0)/k0");

    linalg::SolverConfig solverCfg;
    solverCfg.maxIterations = ctx.linearMaxIterations;
    solverCfg.tolerance = ctx.linearTolerance;
    const auto offlineStarted = std::chrono::steady_clock::now();
    buildOffline(expansionHz, ctx.solver, solverCfg);
    if (options_.fieldPade) {
        const auto lambdaAt = [](double f) {
            const double k = 2.0 * pi * f / c0;
            return k * k;
        };
        const double lambda0 = lambdaAt(expansionHz);
        for (double f : frequencies) {
            const double lambda = lambdaAt(f);
            for (double cutoff : affine_.portCutoffSquared) {
                if ((lambda > cutoff) != (lambda0 > cutoff)) {
                    throw std::runtime_error("AWE field Padé crosses a modal cutoff; split the frequency band");
                }
            }
        }
    }
    const double offlineBuildSec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - offlineStarted).count();

    if (options_.fieldPade) {
        ctx.log.info("AWE paper mode: component-wise field Padé");
    } else {
        ctx.log.info("AWE input Padé order: ["
                     + std::to_string(inputProjectionPade_.numeratorOrder()) + "/"
                     + std::to_string(inputProjectionPade_.denominatorOrder()) + "]");
        ctx.log.info("AWE output Padé order: ["
                     + std::to_string(outputProjectionPade_.numeratorOrder()) + "/"
                     + std::to_string(outputProjectionPade_.denominatorOrder()) + "]");
    }

    SweepResult out;
    out.offlineBuildSec = offlineBuildSec;
    out.points.reserve(frequencies.size());
    const auto onlineStarted = std::chrono::steady_clock::now();
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
    out.onlineSweepSec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - onlineStarted).count();
    if (!ctx.outputDirectory.empty()) {
        fastsweep::FastSweepDiagnostics diag;
        diag.algorithm = name();
        diag.model = options_.fieldPade ? "wavenumber_field_pade" : "wavenumber_port_pade";
        diag.expansionFrequenciesHz = {expansionHz};
        diag.requestedOrder = options_.order;
        if (options_.fieldPade) {
            for (const auto& pade : fieldPade_) {
                diag.padeCoefficientMatchingError = std::max(
                    diag.padeCoefficientMatchingError, pade.coefficientMatchingError());
            }
        } else {
            diag.padeInputPivotRatio = inputProjectionPade_.pivotRatio();
            diag.padeOutputPivotRatio = outputProjectionPade_.pivotRatio();
            diag.padeCoefficientMatchingError = std::max(inputProjectionPade_.coefficientMatchingError(),
                                                          outputProjectionPade_.coefficientMatchingError());
        }
        diag.maxPassivityError = fastsweep::maxPassivityError(out.points);
        const auto path = ctx.outputDirectory / "diagnostics.json";
        if (fastsweep::writeDiagnosticsJson(path, diag)) {
            ctx.log.info("Wrote " + path.string());
        } else {
            ctx.log.warn("Failed to write " + path.string());
        }
    }
    if (options_.fieldPade && !frequencies.empty()) {
        const double t = frequencies.back() / expansionFrequencyHz_ - 1.0;
        out.lastEdgeDofs.reserve(fieldPade_.size());
        for (const auto& pade : fieldPade_) out.lastEdgeDofs.push_back(pade.evaluate(t));
    }
    return out;
}

int AweSweep::buildOffline(double expansionFrequencyHz,
                           linalg::ISparseSolver& solver,
                           const linalg::SolverConfig& solverConfig) {
    if (project_.ports.size() < 2) {
        throw std::runtime_error("AweSweep: at least two project ports are required");
    }

    expansionFrequencyHz_ = expansionFrequencyHz;
    ready_ = false;
    fieldPade_.clear();

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
            momentCount, solver, solverConfig, portVectors,
            fastsweep::PolynomialPortMomentBuilder::Variable::WavenumberShift);

    if (options_.fieldPade) {
        std::vector<Complex> coefficients(fieldMoments.size());
        fieldPade_.reserve(affine_.K.size());
        for (std::size_t i = 0; i < affine_.K.size(); ++i) {
            for (std::size_t r = 0; r < coefficients.size(); ++r) coefficients[r] = fieldMoments[r][i];
            try {
                const bool zero = std::all_of(coefficients.begin(), coefficients.end(),
                    [](Complex value) { return value == Complex{}; });
                fieldPade_.push_back(fastsweep::PadeApproximant::build(coefficients, zero ? 0 : q - 1, zero ? 0 : q));
            } catch (const std::runtime_error& error) {
                throw std::runtime_error("AWE field Padé at DOF " + std::to_string(i) + ": " + error.what());
            }
        }
        ready_ = true;
        return momentCount;
    }

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

    const Complex t(frequencyHz / expansionFrequencyHz_ - 1.0, 0.0);
    if (options_.fieldPade) {
        std::vector<Complex> field;
        field.reserve(fieldPade_.size());
        for (const auto& pade : fieldPade_) field.push_back(pade.evaluate(t));
        return ResultExtractor(project_, portModeSolver_).extract(frequencyHz, field);
    }

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
