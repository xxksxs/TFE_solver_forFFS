#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PolynomialMomentRecurrence.hpp"
#include "bpfem/fastsweep/PortModeUtilities.hpp"
#include "bpfem/linalg/FactorizedSolveSession.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace fem::fastsweep {

namespace {

// 计算广义二项式系数，用于 sqrt(beta) 等非整数幂的局部展开。
double binomialCoefficient(double alpha, int order) {
    double out = 1.0;
    for (int j = 1; j <= order; ++j) {
        out *= (alpha - static_cast<double>(j - 1)) / static_cast<double>(j);
    }
    return out;
}

// 生成 (base + delta)^alpha 关于 delta 的 Taylor 系数。
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

// 对两个截断级数做卷积，保留前 count 阶系数。
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

// 计算小整数幂，避免在矩递推内反复调用通用 pow。
double powInt(double x, int n) {
    double out = 1.0;
    for (int i = 0; i < n; ++i) {
        out *= x;
    }
    return out;
}

}  // namespace

// 将 AffineSystem 中的端口耦合稀疏向量统一提升为全阶复向量。
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

// 在无损材料假设下生成端口边界条件随频率变化的一组 AWE 矩向量。
// 构造 ALPS 文档 P1 方案的一阶端口仿射系统，保留端口导纳与激励的一阶变化。
PolynomialPortMomentBuilder::LosslessLinearization
PolynomialPortMomentBuilder::buildLosslessLinearization(
    const ProjectDefinition& project,
    const FEMAssembler& assembler,
    const PortModeSolver& portModeSolver,
    const FEMAssembler::AffineSystem& affine,
    double expansionFrequencyHz,
    const std::vector<std::vector<Complex>>& portVectors) {
    const int virtualPortCount = static_cast<int>(affine.portCoupling.size());
    if (static_cast<int>(portVectors.size()) != virtualPortCount) {
        throw std::runtime_error(
            "PolynomialPortMomentBuilder: port vector count mismatch in linearization");
    }

    LosslessLinearization out;
    const double kExpansion = 2.0 * pi * expansionFrequencyHz / c0;
    out.lambda0 = kExpansion * kExpansion;
    out.lambdaScale = std::max(out.lambda0, 1.0);
    out.matrixAtExpansion = assembler.assemble(expansionFrequencyHz, out.rhs0);
    out.rhs1.assign(out.rhs0.size(), Complex(0.0, 0.0));
    out.portAdmittanceFirstCoefficients.assign(
        static_cast<std::size_t>(virtualPortCount), 0.0);

    for (int v = 0; v < virtualPortCount; ++v) {
        const double cutoffSquared =
            affine.portCutoffSquared[static_cast<std::size_t>(v)];
        const double betaSquared = out.lambda0 - cutoffSquared;
        if (betaSquared <= 0.0) {
            throw std::runtime_error(
                "Fast-sweep P1 expansion point must be above every retained port-mode cutoff");
        }
        const double beta0 = std::sqrt(betaSquared);
        out.portAdmittanceFirstCoefficients[static_cast<std::size_t>(v)] =
            out.lambdaScale / (2.0 * beta0);

        if (!affine.isExcitationMode[static_cast<std::size_t>(v)]) {
            continue;
        }
        const int projectPort = affine.projectPortIndex[static_cast<std::size_t>(v)];
        if (projectPort < 0 || static_cast<std::size_t>(projectPort) >= project.ports.size()) {
            throw std::runtime_error(
                "PolynomialPortMomentBuilder: invalid project port in linearization");
        }
        const PortMode& mode = virtualPortMode(portModeSolver, affine, v);
        const double normalization0 = powerNormalizationFactor(mode, expansionFrequencyHz);
        if (normalization0 <= 0.0) {
            continue;
        }
        const auto& port = project.ports[static_cast<std::size_t>(projectPort)];
        const Complex incident0 = std::polar(
            std::sqrt(std::max(port.magnitudeW, 0.0)) * normalization0,
            port.phaseDeg * pi / 180.0);
        const Complex rhs0Factor = Complex(0.0, 2.0 * beta0) * incident0;
        const double rhsFirstRatio = 0.25 * out.lambdaScale
            * (1.0 / betaSquared + 1.0 / out.lambda0);
        const Complex rhsFirstFactor = rhs0Factor * rhsFirstRatio;
        const auto& portVector = portVectors[static_cast<std::size_t>(v)];
        for (std::size_t i = 0; i < out.rhs1.size(); ++i) {
            out.rhs1[i] += rhsFirstFactor * portVector[i];
        }
    }
    return out;
}

// 计算 A1*x；A1 同时包含质量矩阵项和端口传播常数的一阶导数项。
std::vector<PolynomialPortMomentBuilder::Complex>
PolynomialPortMomentBuilder::applyLosslessFirstOrderMatrix(
    const FEMAssembler::AffineSystem& affine,
    const LosslessLinearization& linearization,
    const std::vector<std::vector<Complex>>& portVectors,
    const std::vector<Complex>& x) {
    if (x.size() != affine.M.size()
        || portVectors.size() != linearization.portAdmittanceFirstCoefficients.size()) {
        throw std::runtime_error(
            "PolynomialPortMomentBuilder: incompatible vector in first-order matrix action");
    }
    std::vector<Complex> out = affine.M.multiply(x);
    for (auto& value : out) {
        value *= -linearization.lambdaScale;
    }
    for (std::size_t v = 0; v < portVectors.size(); ++v) {
        const double coefficient = linearization.portAdmittanceFirstCoefficients[v];
        const auto& portVector = portVectors[v];
        const Complex projection = bilinear(portVector, x);
        const Complex scale(0.0, coefficient);
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] += scale * projection * portVector[i];
        }
    }
    return out;
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

    linalg::FactorizedSolveSession solveSession(solver, a0, solverConfig);
    auto solveAtExpansion = [&](const std::vector<Complex>& rhs) {
        return solveSession.solve(rhs).field;
    };

    const int virtualPortCount = static_cast<int>(affine.portCoupling.size());
    if (static_cast<int>(portVectors.size()) != virtualPortCount) {
        throw std::runtime_error("PolynomialPortMomentBuilder: port vector count mismatch");
    }

    const std::size_t fullDim = rhs0.size();
    std::vector<std::vector<double>> admittanceCoeffs(
        static_cast<std::size_t>(virtualPortCount),
        std::vector<double>(static_cast<std::size_t>(momentCount), 0.0));
    std::vector<std::vector<Complex>> rhsPortCoefficients(
        static_cast<std::size_t>(virtualPortCount),
        std::vector<Complex>(static_cast<std::size_t>(momentCount), Complex(0.0, 0.0)));

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
            rhsPortCoefficients[static_cast<std::size_t>(v)][static_cast<std::size_t>(r)] =
                coeff;
        }
    }

    // 每一阶只构造当前 RHS，避免保存 momentCount 个全尺寸向量。
    auto rhsCoefficientAt = [&](std::size_t order) {
        if (order == 0) {
            return rhs0;
        }
        std::vector<Complex> rhs(fullDim, Complex(0.0, 0.0));
        for (int v = 0; v < virtualPortCount; ++v) {
            const Complex coeff =
                rhsPortCoefficients[static_cast<std::size_t>(v)][order];
            if (coeff == Complex(0.0, 0.0)) {
                continue;
            }
            const auto& modeVector = portVectors[static_cast<std::size_t>(v)];
            for (std::size_t i = 0; i < fullDim; ++i) {
                rhs[i] += coeff * modeVector[i];
            }
        }
        return rhs;
    };

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
        static_cast<std::size_t>(momentCount), rhsCoefficientAt,
        matrixCoefficientOperators, solveAtExpansion);
}

}  // namespace fem::fastsweep
