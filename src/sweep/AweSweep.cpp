#include "bpfem/sweep/AweSweep.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/PolynomialMomentRecurrence.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem::sweep {

namespace {

using Complex = std::complex<double>;

std::vector<Complex> liftPortVector(const std::vector<std::pair<int, double>>& sparse,
                                    std::size_t n) {
    std::vector<Complex> out(n, Complex(0.0, 0.0));
    for (const auto& [i, w] : sparse) {
        if (i >= 0 && static_cast<std::size_t>(i) < n) {
            out[static_cast<std::size_t>(i)] = w;
        }
    }
    return out;
}

Complex bilinear(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    Complex sum(0.0, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

double binomialCoefficient(double alpha, int order) {
    double out = 1.0;
    for (int j = 1; j <= order; ++j) {
        out *= (alpha - static_cast<double>(j - 1)) / static_cast<double>(j);
    }
    return out;
}

std::vector<double> shiftedPowerSeries(double base, double alpha, int count) {
    if (base <= 0.0) {
        throw std::runtime_error("AweSweep: expansion series base must be positive");
    }
    std::vector<double> coeffs(static_cast<std::size_t>(count), 0.0);
    double invPower = 1.0;
    for (int r = 0; r < count; ++r) {
        coeffs[static_cast<std::size_t>(r)] = binomialCoefficient(alpha, r) * invPower;
        invPower /= base;
    }
    return coeffs;
}

std::vector<double> convolve(const std::vector<double>& a, const std::vector<double>& b, int count) {
    std::vector<double> out(static_cast<std::size_t>(count), 0.0);
    for (int i = 0; i < count; ++i) {
        double sum = 0.0;
        for (int j = 0; j <= i; ++j) {
            sum += a[static_cast<std::size_t>(j)] * b[static_cast<std::size_t>(i - j)];
        }
        out[static_cast<std::size_t>(i)] = sum;
    }
    return out;
}

double powInt(double x, int n) {
    double out = 1.0;
    for (int i = 0; i < n; ++i) {
        out *= x;
    }
    return out;
}

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

    std::vector<Complex> rhs0;
    SparseMatrix A0 = assembler_.assemble(expansionFrequencyHz_, rhs0);
    if (rhs0.empty()) {
        throw std::runtime_error("AweSweep: expansion RHS is empty");
    }

    auto solveAtExpansion = [&](const std::vector<Complex>& rhs) {
        SolveResult result = solver.solve(A0, rhs, solverConfig);
        return result.field;
    };

    const int q = std::max(1, options_.order);
    const int momentCount = 2 * q;
    const std::size_t fullDim = rhs0.size();
    const int virtualPortCount = static_cast<int>(affine_.portCoupling.size());
    std::vector<std::vector<Complex>> portVectors;
    portVectors.reserve(static_cast<std::size_t>(virtualPortCount));
    for (const auto& coupling : affine_.portCoupling) {
        portVectors.push_back(liftPortVector(coupling, fullDim));
    }

    std::vector<std::vector<double>> admittanceCoeffs(
        static_cast<std::size_t>(virtualPortCount),
        std::vector<double>(static_cast<std::size_t>(momentCount), 0.0));
    std::vector<std::vector<Complex>> rhsCoefficients(
        static_cast<std::size_t>(momentCount),
        std::vector<Complex>(fullDim, Complex(0.0, 0.0)));

    for (int v = 0; v < virtualPortCount; ++v) {
        const double kc2 = affine_.portCutoffSquared[static_cast<std::size_t>(v)];
        const double betaBase = lambda0_ - kc2;
        const bool propagatingAtExpansion = betaBase > 0.0;
        const double seriesBase = propagatingAtExpansion ? betaBase : lambda0_;
        const double scalar0 = std::sqrt(seriesBase);
        const auto sqrtSeries = shiftedPowerSeries(seriesBase, 0.5, momentCount);
        for (int r = 0; r < momentCount; ++r) {
            admittanceCoeffs[static_cast<std::size_t>(v)][static_cast<std::size_t>(r)] =
                scalar0 * sqrtSeries[static_cast<std::size_t>(r)] * powInt(lambdaScale_, r);
        }

        if (!affine_.isExcitationMode[static_cast<std::size_t>(v)]) {
            continue;
        }
        const PortMode& mode = virtualPortMode(v);
        const double s0 = powerNormalizationFactor(mode, expansionFrequencyHz_);
        if (s0 <= 0.0) {
            continue;
        }
        const int projectPort = affine_.projectPortIndex[static_cast<std::size_t>(v)];
        const auto& port = project_.ports[static_cast<std::size_t>(projectPort)];
        const Complex incident0 = std::polar(
            std::sqrt(std::max(port.magnitudeW, 0.0)) * s0,
            port.phaseDeg * pi / 180.0);
        const Complex rhsFactor0 = Complex(0.0, 2.0 * admittanceCoeffs[static_cast<std::size_t>(v)][0])
            * incident0;

        // For propagating modes, RHS scalar is proportional to
        // sqrt(beta(lambda) * omega(lambda)).
        const auto betaQuarter = shiftedPowerSeries(seriesBase, 0.25, momentCount);
        const auto lambdaQuarter = shiftedPowerSeries(lambda0_, 0.25, momentCount);
        const auto rhsRatio = propagatingAtExpansion
            ? convolve(betaQuarter, lambdaQuarter, momentCount)
            : std::vector<double>(static_cast<std::size_t>(momentCount), 0.0);
        for (int r = 0; r < momentCount; ++r) {
            const Complex coeff = rhsFactor0
                * rhsRatio[static_cast<std::size_t>(r)]
                * powInt(lambdaScale_, r);
            auto& br = rhsCoefficients[static_cast<std::size_t>(r)];
            const auto& m = portVectors[static_cast<std::size_t>(v)];
            for (std::size_t i = 0; i < fullDim; ++i) {
                br[i] += coeff * m[i];
            }
        }
    }
    rhsCoefficients.front() = rhs0;

    std::vector<fastsweep::PolynomialMomentRecurrence::LinearOperator> matrixCoefficientOperators;
    matrixCoefficientOperators.reserve(static_cast<std::size_t>(momentCount - 1));
    for (int r = 1; r < momentCount; ++r) {
        matrixCoefficientOperators.push_back(
            [this, r, &admittanceCoeffs, &portVectors, virtualPortCount](const std::vector<Complex>& x) {
                std::vector<Complex> out(x.size(), Complex(0.0, 0.0));
                if (r == 1) {
                    std::vector<Complex> mx = affine_.M.multiply(x);
                    for (std::size_t i = 0; i < out.size(); ++i) {
                        out[i] -= lambdaScale_ * mx[i];
                    }
                }
                for (int v = 0; v < virtualPortCount; ++v) {
                    const double coeff =
                        admittanceCoeffs[static_cast<std::size_t>(v)][static_cast<std::size_t>(r)];
                    if (coeff == 0.0) {
                        continue;
                    }
                    const auto& m = portVectors[static_cast<std::size_t>(v)];
                    const Complex projection = bilinear(m, x);
                    const Complex scale(0.0, coeff);
                    for (std::size_t i = 0; i < out.size(); ++i) {
                        out[i] += scale * projection * m[i];
                    }
                }
                return out;
            });
    }

    const auto fieldMoments = fastsweep::PolynomialMomentRecurrence::generatePolynomialMoments(
        rhsCoefficients, matrixCoefficientOperators, solveAtExpansion);

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
    const int faceId = affine_.portFaceIds[static_cast<std::size_t>(virtualPortIndex)];
    if (const auto* multi = portModeSolver_.multiMode(faceId); multi != nullptr) {
        int modeIndex = 0;
        for (int i = 0; i < virtualPortIndex; ++i) {
            if (affine_.portFaceIds[static_cast<std::size_t>(i)] == faceId) {
                ++modeIndex;
            }
        }
        return multi->modes[static_cast<std::size_t>(modeIndex)];
    }
    return portModeSolver_.solve(faceId);
}

}  // namespace fem::sweep
