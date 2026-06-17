#include "bpfem/fastsweep/GalerkinReducedModel.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PortModeUtilities.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace fem::fastsweep {

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

// 在 ROM 空间组装指定频率的小型端口边界系统，并求解 reduced 坐标。
std::vector<GalerkinReducedModel::Complex>
GalerkinReducedModel::solveReduced(double frequencyHz) const {
    if (!ready_ || project_ == nullptr || affine_ == nullptr || portModeSolver_ == nullptr) {
        throw std::runtime_error("GalerkinReducedModel::solveReduced called before build");
    }

    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double k0Sq = k0 * k0;
    std::vector<double> beta(static_cast<std::size_t>(numVirtualPorts_), 0.0);
    std::vector<Complex> incident(static_cast<std::size_t>(numVirtualPorts_), Complex(0.0, 0.0));

    for (int v = 0; v < numVirtualPorts_; ++v) {
        const double kc2 = affine_->portCutoffSquared[static_cast<std::size_t>(v)];
        const double betaSq = k0Sq - kc2;
        beta[static_cast<std::size_t>(v)] = (betaSq > 0.0) ? std::sqrt(betaSq) : 0.0;

        if (affine_->isExcitationMode[static_cast<std::size_t>(v)]) {
            const PortMode& pm = virtualPortMode(*portModeSolver_, *affine_, v);
            const double s = powerNormalizationFactor(pm, frequencyHz);
            if (s > 0.0) {
                const auto& port = project_->ports[static_cast<std::size_t>(
                    affine_->projectPortIndex[static_cast<std::size_t>(v)])];
                incident[static_cast<std::size_t>(v)] = std::polar(
                    std::sqrt(std::max(port.magnitudeW, 0.0)) * s,
                    port.phaseDeg * pi / 180.0);
            }
        }
    }

    const int q = romDim_;
    std::vector<Complex> Atilde(static_cast<std::size_t>(q) * static_cast<std::size_t>(q),
                                Complex(0.0, 0.0));
    for (std::size_t k = 0; k < Atilde.size(); ++k) {
        Atilde[k] = Ktilde_[k] - k0Sq * Mtilde_[k];
    }
    for (int v = 0; v < numVirtualPorts_; ++v) {
        const double bv = beta[static_cast<std::size_t>(v)];
        if (bv <= 0.0) {
            continue;
        }
        const Complex jbeta(0.0, bv);
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
        const double bv = beta[static_cast<std::size_t>(v)];
        if (bv <= 0.0 || incident[static_cast<std::size_t>(v)] == Complex(0.0, 0.0)) {
            continue;
        }
        const Complex factor = Complex(0.0, 2.0) * bv * incident[static_cast<std::size_t>(v)];
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

// 返回当前 ROM 基的最大正交性误差。
double GalerkinReducedModel::basisOrthogonalityError() const {
    return orthogonalityError(basis_);
}

}  // namespace fem::fastsweep
