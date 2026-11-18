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
#include <limits>
#include <vector>

namespace fem::fastsweep {

namespace {

std::vector<double> powerSeries(double constant, double linear, double quadratic,
                                double exponent, int count) {
    std::vector<double> result(static_cast<std::size_t>(count), 0.0);
    result[0] = std::pow(constant, exponent);
    for (int n = 1; n < count; ++n) {
        double coefficient = (exponent + 1.0 - n) * linear * result[n - 1];
        if (n >= 2) {
            coefficient += (2.0 * (exponent + 1.0) - n) * quadratic * result[n - 2];
        }
        result[n] = coefficient / (n * constant);
    }
    return result;
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

double PolynomialPortMomentBuilder::PortTaylorSeries::coordinate(double frequencyHz) const {
    const double ratio = frequencyHz / expansionHz;
    return variable == Variable::WavenumberShift
        ? ratio - 1.0 : (ratio * ratio - 1.0) * lambda0 / scale;
}

PolynomialPortMomentBuilder::PortTaylorSeries PolynomialPortMomentBuilder::buildTaylorSeries(
    const ProjectDefinition& project, const PortModeSolver& modes,
    const FEMAssembler::AffineSystem& affine, double expansionHz, int count,
    Variable variable) {
    if (!std::isfinite(expansionHz) || expansionHz <= 0.0 || count < 1) {
        throw std::invalid_argument("Taylor series requires a positive expansion frequency and count");
    }
    PortTaylorSeries out;
    out.expansionHz = expansionHz;
    out.variable = variable;
    const double k0 = 2.0 * pi * expansionHz / c0;
    out.lambda0 = k0 * k0;
    out.scale = std::max(out.lambda0, 1.0);
    const double linear = variable == Variable::WavenumberShift ? 2.0 * out.lambda0 : out.scale;
    const double quadratic = variable == Variable::WavenumberShift ? out.lambda0 : 0.0;
    out.lambda.assign(static_cast<std::size_t>(count), 0.0);
    out.lambda[0] = out.lambda0;
    if (count > 1) out.lambda[1] = linear;
    if (count > 2) out.lambda[2] = quadratic;
    const auto lambdaQuarter = powerSeries(out.lambda0, linear, quadratic, 0.25, count);

    const auto portCount = affine.portCoupling.size();
    out.admittance.reserve(portCount);
    out.excitation.resize(portCount, std::vector<Complex>(static_cast<std::size_t>(count)));
    for (std::size_t v = 0; v < portCount; ++v) {
        const double betaSquared = out.lambda0 - affine.portCutoffSquared[v];
        if (std::abs(betaSquared) <= 8.0 * std::numeric_limits<double>::epsilon() * out.lambda0) {
            throw std::runtime_error("Cannot Taylor expand a port at its mode cutoff");
        }
        const bool propagating = betaSquared > 0.0;
        const double base = propagating ? betaSquared : out.lambda0;
        out.admittance.push_back(powerSeries(base, linear, quadratic, 0.5, count));
        if (!propagating || !affine.isExcitationMode[v]) continue;
        const auto& port = project.ports.at(static_cast<std::size_t>(affine.projectPortIndex[v]));
        const auto& mode = virtualPortMode(modes, affine, static_cast<int>(v));
        const double norm = powerNormalizationFactor(mode, expansionHz);
        const Complex rhs0 = Complex(0.0, 2.0 * std::sqrt(betaSquared)) *
            std::polar(std::sqrt(std::max(port.magnitudeW, 0.0)) * norm, port.phaseDeg * pi / 180.0);
        const auto betaQuarter = powerSeries(base, linear, quadratic, 0.25, count);
        for (int n = 0; n < count; ++n) {
            double coefficient = 0.0;
            for (int r = 0; r <= n; ++r) {
                coefficient += betaQuarter[r] * lambdaQuarter[n - r];
            }
            out.excitation[v][n] = rhs0 * (coefficient / (betaQuarter[0] * lambdaQuarter[0]));
        }
    }
    return out;
}

std::vector<PolynomialPortMomentBuilder::Complex>
PolynomialPortMomentBuilder::applyMatrixCoefficient(
    const FEMAssembler::AffineSystem& affine, const PortTaylorSeries& series,
    const std::vector<std::vector<Complex>>& ports, std::size_t order,
    const std::vector<Complex>& x) {
    std::vector<Complex> out(x.size());
    const double mass = series.lambda.at(order);
    if (mass != 0.0) {
        out = affine.M.multiply(x);
        for (auto& value : out) value *= -mass;
    }
    if (order == 0) {
        const auto kx = affine.K.multiply(x);
        for (std::size_t i = 0; i < x.size(); ++i) out[i] += kx[i];
        for (int dof : affine.constrainedDofs) out[static_cast<std::size_t>(dof)] = x[static_cast<std::size_t>(dof)];
    }
    for (std::size_t v = 0; v < ports.size(); ++v) {
        const auto& entries = affine.portCoupling[v];
        const Complex coefficient(0.0, series.admittance[v][order]);
        Complex projection = 0.0;
        for (const auto& [dof, weight] : entries) projection += weight * x[static_cast<std::size_t>(dof)];
        for (const auto& [dof, weight] : entries) out[static_cast<std::size_t>(dof)] += coefficient * projection * weight;
    }
    return out;
}

std::vector<PolynomialPortMomentBuilder::Complex>
PolynomialPortMomentBuilder::rhsCoefficient(
    const PortTaylorSeries& series, const std::vector<std::vector<Complex>>& ports,
    std::size_t dimension, std::size_t order) {
    std::vector<Complex> out(dimension);
    for (std::size_t v = 0; v < ports.size(); ++v) {
        const auto coefficient = series.excitation[v].at(order);
        if (coefficient == Complex{}) continue;
        for (std::size_t i = 0; i < dimension; ++i) out[i] += coefficient * ports[v][i];
    }
    return out;
}

std::vector<std::vector<PolynomialPortMomentBuilder::Complex>>
PolynomialPortMomentBuilder::generateLosslessMoments(
    const ProjectDefinition& project, const FEMAssembler& assembler,
    const PortModeSolver& portModeSolver, const FEMAssembler::AffineSystem& affine,
    double expansionFrequencyHz, int momentCount, linalg::ISparseSolver& solver,
    const linalg::SolverConfig& config, const std::vector<std::vector<Complex>>& ports,
    Variable variable,
    const std::function<bool(const SparseMatrix&, const std::vector<Complex>&)>& accept,
    const std::function<void(const std::vector<Complex>&)>& onMoment) {
    const auto series = buildTaylorSeries(project, portModeSolver, affine, expansionFrequencyHz,
                                         momentCount, variable);
    std::vector<Complex> rhs0;
    const auto a0 = assembler.assemble(expansionFrequencyHz, rhs0);
    linalg::FactorizedSolveSession session(solver, a0, config);
    std::size_t solvedMoments = 0;
    std::vector<PolynomialMomentRecurrence::LinearOperator> operators;
    for (int r = 1; r < momentCount; ++r) {
        operators.push_back([&, r](const std::vector<Complex>& x) {
            return applyMatrixCoefficient(affine, series, ports, static_cast<std::size_t>(r), x);
        });
    }
    return PolynomialMomentRecurrence::generatePolynomialMoments(
        static_cast<std::size_t>(momentCount),
        [&](std::size_t r) { return r == 0 ? rhs0 : rhsCoefficient(series, ports, rhs0.size(), r); },
        operators, [&](const std::vector<Complex>& rhs) {
            auto result = session.solve(rhs);
            if (result.status != SolveStatus::Converged) {
                const auto context = std::string(solver.name()) + " at expansion "
                    + std::to_string(expansionFrequencyHz) + " Hz, moment "
                    + std::to_string(solvedMoments);
                linalg::requireConverged(result, context.c_str());
            }
            ++solvedMoments;
            if (onMoment) onMoment(result.field);
            return std::move(result.field);
        }, [&](const std::vector<Complex>& rhs) { return !accept || accept(a0, rhs); });
}

}  // namespace fem::fastsweep
