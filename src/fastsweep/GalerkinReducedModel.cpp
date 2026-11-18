#include "bpfem/fastsweep/GalerkinReducedModel.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/bc/WavePortBC.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PortModeUtilities.hpp"

#include <algorithm>
#include <limits>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <algorithm>
#include <limits>

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
}  // namespace

void GalerkinReducedModel::setPolynomials(
    std::vector<PolynomialPortMomentBuilder::PortTaylorSeries> polynomials, bool lastExpansionOnly) {
    polynomials_ = std::move(polynomials);
    lastExpansionOnly_ = lastExpansionOnly;
}

// 用给定正交基构造 Galerkin ROM：投影 K/M 矩阵并预投影端口向量。
void GalerkinReducedModel::build(
    const ProjectDefinition& project,
    const FEMAssembler::AffineSystem& affine,
    const PortModeSolver& portModeSolver,
    std::vector<std::vector<Complex>> basis,
    const std::vector<std::vector<Complex>>& portVectors) {
    if (basis.empty()) {
        throw std::runtime_error("GalerkinReducedModel: basis is empty");
    }
    polynomials_.clear();
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
        const auto Kvj = affine.K.multiply(basis_[static_cast<std::size_t>(j)]);
        const auto Mvj = affine.M.multiply(basis_[static_cast<std::size_t>(j)]);
        for (int i = 0; i < romDim_; ++i) {
            Ktilde_[static_cast<std::size_t>(i) * static_cast<std::size_t>(romDim_)
                    + static_cast<std::size_t>(j)] =
                bilinear(basis_[static_cast<std::size_t>(i)], Kvj);
            Mtilde_[static_cast<std::size_t>(i) * static_cast<std::size_t>(romDim_)
                    + static_cast<std::size_t>(j)] =
                bilinear(basis_[static_cast<std::size_t>(i)], Mvj);
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

const PolynomialPortMomentBuilder::PortTaylorSeries*
GalerkinReducedModel::polynomialAt(double frequencyHz) const {
    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double k0Sq = k0 * k0;
    const PolynomialPortMomentBuilder::PortTaylorSeries* series = nullptr;
    for (const auto& candidate : polynomials_) {
        if (lastExpansionOnly_ && &candidate != &polynomials_.back()) continue;
        const bool sameBranch = std::all_of(affine_->portCutoffSquared.begin(),
            affine_->portCutoffSquared.end(), [&](double cutoff) {
                return (k0Sq > cutoff) == (candidate.lambda0 > cutoff);
            });
        if (sameBranch && (!series || std::abs(frequencyHz - candidate.expansionHz)
                                   < std::abs(frequencyHz - series->expansionHz))) {
            series = &candidate;
        }
    }
    return series;
}

GalerkinReducedModel::Coefficients GalerkinReducedModel::coefficientsAt(double frequencyHz) const {
    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double k0Sq = k0 * k0;
    const auto* series = polynomialAt(frequencyHz);
    if (!polynomials_.empty() && !series) {
        throw std::runtime_error("No polynomial expansion on the requested port cutoff branch");
    }
    const double coordinate = series ? series->coordinate(frequencyHz) : 0.0;
    Coefficients coefficients;
    coefficients.mass = series ? evaluateSeries(series->lambda, coordinate) : k0Sq;
    coefficients.admittance.resize(static_cast<std::size_t>(numVirtualPorts_));
    coefficients.excitation.resize(static_cast<std::size_t>(numVirtualPorts_));
    for (int v = 0; v < numVirtualPorts_; ++v) {
        const auto index = static_cast<std::size_t>(v);
        if (series) {
            coefficients.admittance[index] = evaluateSeries(series->admittance[index], coordinate);
            coefficients.excitation[index] = evaluateSeries(series->excitation[index], coordinate);
            continue;
        }
        const double cutoff = affine_->portCutoffSquared[index];
        const double beta = bc::wavePortAdmittance(k0Sq, cutoff);
        coefficients.admittance[index] = beta;
        if (k0Sq > cutoff && affine_->isExcitationMode[index]) {
            const auto& mode = virtualPortMode(*portModeSolver_, *affine_, v);
            const auto& port = project_->ports[static_cast<std::size_t>(affine_->projectPortIndex[index])];
            const double amplitude = std::sqrt(std::max(port.magnitudeW, 0.0))
                * powerNormalizationFactor(mode, frequencyHz);
            coefficients.excitation[index] = Complex(0.0, 2.0 * beta)
                * std::polar(amplitude, port.phaseDeg * pi / 180.0);
        }
    }
    return coefficients;
}

// 在 ROM 空间组装指定频率的小型端口边界系统，并求解 reduced 坐标。
std::vector<GalerkinReducedModel::Complex>
GalerkinReducedModel::solveReduced(double frequencyHz) const {
    if (!ready_ || project_ == nullptr || affine_ == nullptr || portModeSolver_ == nullptr) {
        throw std::runtime_error("GalerkinReducedModel::solveReduced called before build");
    }

    const auto coefficients = coefficientsAt(frequencyHz);
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
            const Complex jbVi = jbeta * vp[static_cast<std::size_t>(i)];
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
            rhs[static_cast<std::size_t>(i)] += factor * vp[static_cast<std::size_t>(i)];
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
SParameterPoint GalerkinReducedModel::evaluate(double frequencyHz) const {
    SParameterPoint sp;
    sp.frequencyHz = frequencyHz;
    if (!ready_ || numProjectPorts_ < 2 || affine_ == nullptr || project_ == nullptr
        || portModeSolver_ == nullptr) {
        return sp;
    }

    const auto xTilde = solveReduced(frequencyHz);
    std::vector<double> sNorm(static_cast<std::size_t>(numVirtualPorts_), 0.0);
    for (int v = 0; v < numVirtualPorts_; ++v) {
        const PortMode& pm = virtualPortMode(*portModeSolver_, *affine_, v);
        sNorm[static_cast<std::size_t>(v)] = powerNormalizationFactor(pm, frequencyHz);
    }
    const auto dominant = dominantVirtualPortsByProject(*portModeSolver_, *affine_,
                                                        numProjectPorts_);

    auto portProjection = [&](int virtualIdx) {
        const auto& vp = portModeReduced_[static_cast<std::size_t>(virtualIdx)];
        Complex sum(0.0, 0.0);
        for (int i = 0; i < romDim_; ++i) {
            sum += vp[static_cast<std::size_t>(i)] * xTilde[static_cast<std::size_t>(i)];
        }
        return sum;
    };

    const int vIn = dominant[0];
    const int vOut = (numProjectPorts_ >= 2) ? dominant[1] : -1;
    if (vIn < 0) {
        return sp;
    }
    const double sIn = sNorm[static_cast<std::size_t>(vIn)];
    const double sOut = (vOut >= 0) ? sNorm[static_cast<std::size_t>(vOut)] : 0.0;
    if (sIn <= 0.0) {
        return sp;
    }

    const Complex bIn = portProjection(vIn) / sIn;
    const Complex bOut = (vOut >= 0 && sOut > 0.0)
        ? portProjection(vOut) / sOut
        : Complex(0.0, 0.0);
    const auto& port0 = project_->ports[0];
    const Complex incidentExtract = port0.excited
        ? std::polar(std::sqrt(std::max(port0.magnitudeW, 0.0)),
                     port0.phaseDeg * pi / 180.0)
        : Complex(1.0, 0.0);
    if (incidentExtract == Complex(0.0, 0.0)) {
        return sp;
    }
    sp.s11 = (bIn - incidentExtract) / incidentExtract;
    sp.s21 = bOut / incidentExtract;
    return sp;
}

// 将 reduced 坐标按基向量线性组合回全阶场，用于 VTU 场文件输出。
std::vector<GalerkinReducedModel::Complex>
GalerkinReducedModel::reconstructField(double frequencyHz) const {
    const auto xTilde = solveReduced(frequencyHz);
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

double GalerkinReducedModel::relativeResidual(double frequencyHz) const {
    if (!polynomials_.empty() && !polynomialAt(frequencyHz)) return std::numeric_limits<double>::infinity();
    const auto x = reconstructField(frequencyHz);
    const auto coefficients = coefficientsAt(frequencyHz);
    auto residual = affine_->K.multiply(x);
    const auto mx = affine_->M.multiply(x);
    for (std::size_t i = 0; i < x.size(); ++i) residual[i] -= coefficients.mass * mx[i];
    std::vector<Complex> rhs(x.size());
    for (int dof : affine_->constrainedDofs) residual[static_cast<std::size_t>(dof)] = x[static_cast<std::size_t>(dof)];
    for (std::size_t v = 0; v < affine_->portCoupling.size(); ++v) {
        Complex projection{};
        for (const auto& [dof, weight] : affine_->portCoupling[v]) {
            projection += weight * x[static_cast<std::size_t>(dof)];
        }
        const Complex scale = Complex(0.0, coefficients.admittance[v]) * projection;
        for (const auto& [dof, weight] : affine_->portCoupling[v]) {
            residual[static_cast<std::size_t>(dof)] += scale * weight;
            rhs[static_cast<std::size_t>(dof)] += coefficients.excitation[v] * weight;
        }
    }
    double residualMax = 0.0;
    double rhsMax = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        residualMax = std::max(residualMax, std::abs(residual[i] - rhs[i]));
        rhsMax = std::max(rhsMax, std::abs(rhs[i]));
    }
    return residualMax / std::max(rhsMax, 1.0e-300);
}

// 返回当前 ROM 基的最大正交性误差。
double GalerkinReducedModel::basisOrthogonalityError() const {
    return orthogonalityError(basis_);
}

}  // namespace fem::fastsweep
