#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PolynomialMomentRecurrence.hpp"
#include "bpfem/fastsweep/PortModeUtilities.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace fem::fastsweep {

namespace {

double binomialCoefficient(double alpha, int order) {
    double out = 1.0;
    for (int j = 1; j <= order; ++j) {
        out *= (alpha - static_cast<double>(j - 1)) / static_cast<double>(j);
    }
    return out;
}

std::vector<double> shiftedPowerSeries(double base, double alpha, int count) {
    if (base <= 0.0) {
        throw std::runtime_error("PolynomialPortMomentBuilder: expansion series base must be positive");
    }
    std::vector<double> coeffs(static_cast<std::size_t>(count), 0.0);
    double invPower = 1.0;
    for (int r = 0; r < count; ++r) {
        coeffs[static_cast<std::size_t>(r)] = binomialCoefficient(alpha, r) * invPower;
        invPower /= base;
    }
    return coeffs;
}

std::vector<double> convolve(const std::vector<double>& a,
                             const std::vector<double>& b,
                             int count) {
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

std::vector<std::vector<PolynomialPortMomentBuilder::Complex>>
PolynomialPortMomentBuilder::buildPortVectors(const FEMAssembler::AffineSystem& affine,
                                              std::size_t fullDimension) {
    std::vector<std::vector<Complex>> portVectors;
    portVectors.reserve(affine.portCoupling.size());
    for (const auto& coupling : affine.portCoupling) {
        portVectors.push_back(liftRealSparseVector(coupling, fullDimension));
    }
    return portVectors;
}

std::vector<std::vector<PolynomialPortMomentBuilder::Complex>>
PolynomialPortMomentBuilder::generateLosslessMoments(
    const ProjectDefinition& project,
    const FEMAssembler& assembler,
    const PortModeSolver& portModeSolver,
    const FEMAssembler::AffineSystem& affine,
    double expansionFrequencyHz,
    int momentCount,
    linalg::ISparseSolver& solver,
    const linalg::SolverConfig& solverConfig,
    const std::vector<std::vector<Complex>>& portVectors) {
    if (momentCount < 1) {
        throw std::runtime_error("PolynomialPortMomentBuilder: momentCount must be >= 1");
    }

    const double kExpansion = 2.0 * pi * expansionFrequencyHz / c0;
    const double lambda0 = kExpansion * kExpansion;
    const double lambdaScale = std::max(lambda0, 1.0);

    std::vector<Complex> rhs0;
    SparseMatrix a0 = assembler.assemble(expansionFrequencyHz, rhs0);
    if (rhs0.empty()) {
        throw std::runtime_error("PolynomialPortMomentBuilder: expansion RHS is empty");
    }

    auto solveAtExpansion = [&](const std::vector<Complex>& rhs) {
        SolveResult result = solver.solve(a0, rhs, solverConfig);
        return result.field;
    };

    const int virtualPortCount = static_cast<int>(affine.portCoupling.size());
    if (static_cast<int>(portVectors.size()) != virtualPortCount) {
        throw std::runtime_error("PolynomialPortMomentBuilder: port vector count mismatch");
    }

    const std::size_t fullDim = rhs0.size();
    std::vector<std::vector<double>> admittanceCoeffs(
        static_cast<std::size_t>(virtualPortCount),
        std::vector<double>(static_cast<std::size_t>(momentCount), 0.0));
    std::vector<std::vector<Complex>> rhsCoefficients(
        static_cast<std::size_t>(momentCount),
        std::vector<Complex>(fullDim, Complex(0.0, 0.0)));

    for (int v = 0; v < virtualPortCount; ++v) {
        const double kc2 = affine.portCutoffSquared[static_cast<std::size_t>(v)];
        const double betaBase = lambda0 - kc2;
        const bool propagatingAtExpansion = betaBase > 0.0;
        const double seriesBase = propagatingAtExpansion ? betaBase : lambda0;
        const double scalar0 = std::sqrt(seriesBase);
        const auto sqrtSeries = shiftedPowerSeries(seriesBase, 0.5, momentCount);
        for (int r = 0; r < momentCount; ++r) {
            admittanceCoeffs[static_cast<std::size_t>(v)][static_cast<std::size_t>(r)] =
                scalar0 * sqrtSeries[static_cast<std::size_t>(r)] * powInt(lambdaScale, r);
        }

        if (!affine.isExcitationMode[static_cast<std::size_t>(v)]) {
            continue;
        }

        const PortMode& mode = virtualPortMode(portModeSolver, affine, v);
        const double s0 = powerNormalizationFactor(mode, expansionFrequencyHz);
        if (s0 <= 0.0) {
            continue;
        }
        const int projectPort = affine.projectPortIndex[static_cast<std::size_t>(v)];
        const auto& port = project.ports[static_cast<std::size_t>(projectPort)];
        const Complex incident0 = std::polar(
            std::sqrt(std::max(port.magnitudeW, 0.0)) * s0,
            port.phaseDeg * pi / 180.0);
        const Complex rhsFactor0 =
            Complex(0.0, 2.0 * admittanceCoeffs[static_cast<std::size_t>(v)][0])
            * incident0;

        const auto betaQuarter = shiftedPowerSeries(seriesBase, 0.25, momentCount);
        const auto lambdaQuarter = shiftedPowerSeries(lambda0, 0.25, momentCount);
        const auto rhsRatio = propagatingAtExpansion
            ? convolve(betaQuarter, lambdaQuarter, momentCount)
            : std::vector<double>(static_cast<std::size_t>(momentCount), 0.0);
        for (int r = 0; r < momentCount; ++r) {
            const Complex coeff = rhsFactor0
                * rhsRatio[static_cast<std::size_t>(r)]
                * powInt(lambdaScale, r);
            auto& br = rhsCoefficients[static_cast<std::size_t>(r)];
            const auto& m = portVectors[static_cast<std::size_t>(v)];
            for (std::size_t i = 0; i < fullDim; ++i) {
                br[i] += coeff * m[i];
            }
        }
    }
    rhsCoefficients.front() = rhs0;

    std::vector<PolynomialMomentRecurrence::LinearOperator> matrixCoefficientOperators;
    matrixCoefficientOperators.reserve(static_cast<std::size_t>(std::max(0, momentCount - 1)));
    for (int r = 1; r < momentCount; ++r) {
        matrixCoefficientOperators.push_back(
            [r, lambdaScale, &affine, &admittanceCoeffs, &portVectors, virtualPortCount](
                const std::vector<Complex>& x) {
                std::vector<Complex> out(x.size(), Complex(0.0, 0.0));
                if (r == 1) {
                    std::vector<Complex> mx = affine.M.multiply(x);
                    for (std::size_t i = 0; i < out.size(); ++i) {
                        out[i] -= lambdaScale * mx[i];
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

    return PolynomialMomentRecurrence::generatePolynomialMoments(
        rhsCoefficients, matrixCoefficientOperators, solveAtExpansion);
}

}  // namespace fem::fastsweep
