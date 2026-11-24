#include "bpfem/fastsweep/GalerkinReducedModel.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/bc/WavePortBC.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PortModeUtilities.hpp"
#include "bpfem/post/ResultExtractor.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace fem::fastsweep {

namespace {
template<class Scalar>
Scalar evaluateSeries(const std::vector<Scalar>& coefficients, double coordinate) {
    Scalar result{};
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
        result = result * coordinate + *it;
    }
    return result;
}

using Complex = std::complex<double>;
using PortTaylorSeries = PolynomialPortMomentBuilder::PortTaylorSeries;

struct Coefficients {
    double mass = 0.0;
    std::vector<double> admittance;
    std::vector<Complex> excitation;
};


Coefficients coefficientsAt(const ProjectDefinition& project, const PortModeSolver& modes,
    const FEMAssembler::AffineSystem& affine, double frequencyHz, const PortTaylorSeries* series) {
    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double k0Sq = k0 * k0;
    const double coordinate = series ? series->coordinate(frequencyHz) : 0.0;
    Coefficients coefficients;
    coefficients.mass = series ? evaluateSeries(series->lambda, coordinate) : k0Sq;
    const auto portCount = affine.portCoupling.size();
    coefficients.admittance.resize(portCount);
    coefficients.excitation.resize(portCount);
    for (std::size_t v = 0; v < portCount; ++v) {
        const auto index = static_cast<std::size_t>(v);
        if (series) {
            coefficients.admittance[index] = evaluateSeries(series->admittance[index], coordinate);
            coefficients.excitation[index] = evaluateSeries(series->excitation[index], coordinate);
            continue;
        }
        const double cutoff = affine.portCutoffSquared[index];
        const double beta = bc::wavePortAdmittance(k0Sq, cutoff);
        coefficients.admittance[index] = beta;
        if (k0Sq > cutoff && affine.isExcitationMode[index]) {
            const auto& mode = virtualPortMode(modes, affine, static_cast<int>(v));
            const auto& port = project.ports[static_cast<std::size_t>(affine.projectPortIndex[index])];
            const double amplitude = std::sqrt(std::max(port.magnitudeW, 0.0))
                * powerNormalizationFactor(mode, frequencyHz);
            coefficients.excitation[index] = Complex(0.0, 2.0 * beta)
                * std::polar(amplitude, port.phaseDeg * pi / 180.0);
        }
    }
    return coefficients;
}

struct ResidualNorms {
    double absolute = 0.0;
    double rhs = 0.0;
};

std::vector<Complex> assembleRhs(const FEMAssembler::AffineSystem& affine,
                                 const Coefficients& coefficients) {
    std::vector<Complex> rhs(affine.K.size());
    for (std::size_t port = 0; port < affine.portCoupling.size(); ++port) {
        if (coefficients.excitation[port] == Complex{}) continue;
        for (const auto& [dof, weight] : affine.portCoupling[port]) {
            rhs[static_cast<std::size_t>(dof)] += coefficients.excitation[port] * weight;
        }
    }
    return rhs;
}

ResidualNorms residualNorms(const FEMAssembler::AffineSystem& affine,
                            const std::vector<Complex>& field,
                            const std::vector<Complex>& stiffness,
                            const std::vector<Complex>& mass,
                            const Coefficients& coefficients) {
    auto applied = stiffness;
    for (std::size_t i = 0; i < field.size(); ++i) applied[i] -= coefficients.mass * mass[i];
    for (int dof : affine.constrainedDofs) {
        const auto index = static_cast<std::size_t>(dof);
        applied[index] = field[index];
    }
    for (std::size_t port = 0; port < affine.portCoupling.size(); ++port) {
        Complex projection{};
        for (const auto& [dof, weight] : affine.portCoupling[port]) {
            projection += weight * field[static_cast<std::size_t>(dof)];
        }
        const Complex scale = Complex(0.0, coefficients.admittance[port]) * projection;
        for (const auto& [dof, weight] : affine.portCoupling[port]) {
            applied[static_cast<std::size_t>(dof)] += scale * weight;
        }
    }
    const auto rhs = assembleRhs(affine, coefficients);
    ResidualNorms norms;
    for (std::size_t i = 0; i < field.size(); ++i) {
        const double defect = std::abs(rhs[i] - applied[i]);
        const double amplitude = std::abs(rhs[i]);
        if (!std::isfinite(defect) || !std::isfinite(amplitude)) {
            return {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        }
        norms.absolute = std::max(norms.absolute, defect);
        norms.rhs = std::max(norms.rhs, amplitude);
    }
    return norms;
}

}  // namespace

bool nearPortCutoff(const FEMAssembler::AffineSystem& affine, double frequencyHz) {
    const double k = 2.0 * pi * frequencyHz / c0;
    const double lambda = k * k;
    constexpr double relativeMargin = 1.0e-8;
    return std::any_of(affine.portCutoffSquared.begin(), affine.portCutoffSquared.end(),
        [lambda](double cutoff) {
            return std::abs(lambda - cutoff) <= relativeMargin * std::max(lambda, cutoff);
        });
}

bool hasLosslessExcitation(const ProjectDefinition& project, const PortModeSolver& modes,
                           const FEMAssembler::AffineSystem& affine, double frequencyHz) {
    const auto rhs = assembleRhs(affine, coefficientsAt(project, modes, affine, frequencyHz, nullptr));
    return std::any_of(rhs.begin(), rhs.end(), [](Complex value) { return value != Complex{}; });
}

SweepResidualSample losslessResiduals(
    const ProjectDefinition& project, const PortModeSolver& modes,
    const FEMAssembler::AffineSystem& affine, double frequencyHz,
    const std::vector<Complex>& field, const PortTaylorSeries* series) {
    SweepResidualSample sample;
    sample.frequencyHz = frequencyHz;
    sample.nearCutoff = nearPortCutoff(affine, frequencyHz);
    const auto stiffness = affine.K.multiply(field);
    const auto mass = affine.M.multiply(field);
    const auto exact = coefficientsAt(project, modes, affine, frequencyHz, nullptr);
    const auto trueNorms = residualNorms(affine, field, stiffness, mass, exact);
    sample.zeroRhs = trueNorms.rhs == 0.0;
    if (std::isfinite(trueNorms.absolute)) {
        sample.trueAbsolute = trueNorms.absolute;
        if (!sample.zeroRhs) sample.trueRelative = trueNorms.absolute / trueNorms.rhs;
    }
    const bool sameBranch = series && std::all_of(affine.portCutoffSquared.begin(),
        affine.portCutoffSquared.end(), [&](double cutoff) {
            return (exact.mass > cutoff) == (series->lambda0 > cutoff);
        });
    if (sameBranch && !sample.nearCutoff) {
        const auto taylor = coefficientsAt(project, modes, affine, frequencyHz, series);
        const auto norms = residualNorms(affine, field, stiffness, mass, taylor);
        if (norms.rhs > 0.0 && std::isfinite(norms.absolute)) {
            sample.taylorRelative = norms.absolute / norms.rhs;
        }
    }
    return sample;
}

// 用给定正交基构造 Galerkin ROM：投影 K/M 矩阵并预投影端口向量。
void GalerkinReducedModel::build(
    const ProjectDefinition& project,
    const FEMAssembler::AffineSystem& affine,
    const PortModeSolver& portModeSolver,
    std::vector<std::vector<Complex>> basis,
    const std::vector<std::vector<Complex>>& portVectors,
    Projection projection) {
    if (basis.empty()) {
        throw std::runtime_error("GalerkinReducedModel: basis is empty");
    }
    lastSolveSucceeded_ = true;
    projection_ = projection;
    const auto testProduct = projection == Projection::Hermitian ? hdot : bilinear;
    project_ = &project;
    affine_ = &affine;
    portModeSolver_ = &portModeSolver;
    basis_ = std::move(basis);

    romDim_ = static_cast<int>(basis_.size());
    fullDim_ = static_cast<int>(affine.K.size());
    numProjectPorts_ = static_cast<int>(project.ports.size());
    numVirtualPorts_ = static_cast<int>(affine.portCoupling.size());
    if (static_cast<int>(portVectors.size()) != numVirtualPorts_) {
        throw std::runtime_error("GalerkinReducedModel: port vector count mismatch");
    }

    Ktilde_.assign(static_cast<std::size_t>(romDim_) * static_cast<std::size_t>(romDim_),
                   Complex(0.0, 0.0));
    Mtilde_.assign(static_cast<std::size_t>(romDim_) * static_cast<std::size_t>(romDim_),
                   Complex(0.0, 0.0));

    for (int j = 0; j < romDim_; ++j) {
        const auto& column = basis_[static_cast<std::size_t>(j)];
        auto Kvj = affine.K.multiply(column);
        for (int dof : affine.constrainedDofs) {
            const auto index = static_cast<std::size_t>(dof);
            Kvj[index] = column[index];
        }
        const auto Mvj = affine.M.multiply(column);
        for (int i = 0; i < romDim_; ++i) {
            Ktilde_[static_cast<std::size_t>(i) * static_cast<std::size_t>(romDim_)
                    + static_cast<std::size_t>(j)] =
                testProduct(basis_[static_cast<std::size_t>(i)], Kvj);
            Mtilde_[static_cast<std::size_t>(i) * static_cast<std::size_t>(romDim_)
                    + static_cast<std::size_t>(j)] =
                testProduct(basis_[static_cast<std::size_t>(i)], Mvj);
        }
    }

    portModeReduced_.assign(static_cast<std::size_t>(numVirtualPorts_), {});
    for (int v = 0; v < numVirtualPorts_; ++v) {
        portModeReduced_[static_cast<std::size_t>(v)].assign(static_cast<std::size_t>(romDim_),
                                                              Complex(0.0, 0.0));
        for (int i = 0; i < romDim_; ++i) {
            portModeReduced_[static_cast<std::size_t>(v)][static_cast<std::size_t>(i)] =
                bilinear(basis_[static_cast<std::size_t>(i)],
                         portVectors[static_cast<std::size_t>(v)]);
        }
    }

    ready_ = true;
}

// 在 ROM 空间组装指定频率的小型端口边界系统，并求解 reduced 坐标。
std::vector<GalerkinReducedModel::Complex>
GalerkinReducedModel::solveReduced(double frequencyHz, const PortTaylorSeries* series) const {
    if (!ready_ || project_ == nullptr || affine_ == nullptr || portModeSolver_ == nullptr) {
        throw std::runtime_error("GalerkinReducedModel::solveReduced called before build");
    }

    const auto coefficients = coefficientsAt(*project_, *portModeSolver_, *affine_, frequencyHz, series);
    const int q = romDim_;
    std::vector<Complex> Atilde(static_cast<std::size_t>(q) * static_cast<std::size_t>(q),
                                Complex(0.0, 0.0));
    for (std::size_t k = 0; k < Atilde.size(); ++k) {
        Atilde[k] = Ktilde_[k] - coefficients.mass * Mtilde_[k];
    }
    for (int v = 0; v < numVirtualPorts_; ++v) {
        const Complex jbeta(0.0, coefficients.admittance[static_cast<std::size_t>(v)]);
        const auto& vp = portModeReduced_[static_cast<std::size_t>(v)];
        for (int i = 0; i < q; ++i) {
            const Complex right = vp[static_cast<std::size_t>(i)];
            // The coupling vectors are real, so V^H p = conj(p^T V).
            const Complex left = projection_ == Projection::Hermitian ? std::conj(right) : right;
            const Complex jbVi = jbeta * left;
            for (int j = 0; j < q; ++j) {
                Atilde[static_cast<std::size_t>(i) * static_cast<std::size_t>(q)
                       + static_cast<std::size_t>(j)] +=
                    jbVi * vp[static_cast<std::size_t>(j)];
            }
        }
    }

    std::vector<Complex> rhs(static_cast<std::size_t>(q), Complex(0.0, 0.0));
    for (int v = 0; v < numVirtualPorts_; ++v) {
        const Complex factor = coefficients.excitation[static_cast<std::size_t>(v)];
        const auto& vp = portModeReduced_[static_cast<std::size_t>(v)];
        for (int i = 0; i < q; ++i) {
            const Complex right = vp[static_cast<std::size_t>(i)];
            const Complex left = projection_ == Projection::Hermitian ? std::conj(right) : right;
            rhs[static_cast<std::size_t>(i)] += factor * left;
        }
    }

    std::vector<Complex> xTilde = rhs;
    if (!denseSolve(Atilde, xTilde, q)) {
        lastSolveSucceeded_ = false;
        throw std::runtime_error("GalerkinReducedModel::solveReduced: reduced system singular at "
                                 + std::to_string(frequencyHz) + " Hz");
    }
    lastSolveSucceeded_ = true;
    return xTilde;
}

// 在指定频率求解 ROM 后，从主导输入/输出端口投影中提取 S11 和 S21。
SParameterPoint GalerkinReducedModel::evaluate(double frequencyHz, const PortTaylorSeries* series) const {
    SParameterPoint sp;
    sp.frequencyHz = frequencyHz;
    if (!ready_ || numProjectPorts_ < 2 || affine_ == nullptr || project_ == nullptr
        || portModeSolver_ == nullptr) {
        return sp;
    }

    const auto xTilde = solveReduced(frequencyHz, series);
    const auto dominant = dominantVirtualPortsByProject(*portModeSolver_, *affine_,
                                                        numProjectPorts_);

    auto portProjection = [&](int virtualIdx) {
        if (virtualIdx < 0) return Complex{};
        const auto& vp = portModeReduced_[static_cast<std::size_t>(virtualIdx)];
        return bilinear(vp, xTilde);
    };

    const int vIn = dominant[0];
    const int vOut = (numProjectPorts_ >= 2) ? dominant[1] : -1;
    if (vIn < 0) {
        return sp;
    }
    return ResultExtractor(*project_, *portModeSolver_).extractFromPortProjections(
        frequencyHz, portProjection(vIn), portProjection(vOut));
}

// 将 reduced 坐标按基向量线性组合回全阶场，用于 VTU 场文件输出。
std::vector<GalerkinReducedModel::Complex>
GalerkinReducedModel::reconstructField(double frequencyHz, const PortTaylorSeries* series) const {
    const auto xTilde = solveReduced(frequencyHz, series);
    std::vector<Complex> field(static_cast<std::size_t>(fullDim_), Complex(0.0, 0.0));
    for (int j = 0; j < romDim_; ++j) {
        const auto& v = basis_[static_cast<std::size_t>(j)];
        const Complex coeff = xTilde[static_cast<std::size_t>(j)];
        for (int i = 0; i < fullDim_; ++i) {
            field[static_cast<std::size_t>(i)] += coeff * v[static_cast<std::size_t>(i)];
        }
    }
    return field;
}

double GalerkinReducedModel::relativeResidual(double frequencyHz, const PortTaylorSeries* series) const {
    const auto x = reconstructField(frequencyHz, series);
    const auto coefficients = coefficientsAt(*project_, *portModeSolver_, *affine_, frequencyHz, series);
    const auto norms = residualNorms(*affine_, x, affine_->K.multiply(x),
                                     affine_->M.multiply(x), coefficients);
    return norms.absolute / std::max(norms.rhs, 1.0e-300);
}

// 返回当前 ROM 基的最大正交性误差。
double GalerkinReducedModel::basisOrthogonalityError() const {
    return orthogonalityError(basis_);
}

}  // namespace fem::fastsweep
