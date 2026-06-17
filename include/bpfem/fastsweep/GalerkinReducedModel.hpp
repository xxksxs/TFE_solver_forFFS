#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"

#include <complex>
#include <vector>

namespace fem::fastsweep {

class GalerkinReducedModel {
public:
    using Complex = std::complex<double>;

    // 构造 Galerkin ROM，并缓存投影后的 K/M 矩阵和端口向量。
    void build(const ProjectDefinition& project,
               const FEMAssembler::AffineSystem& affine,
               const PortModeSolver& portModeSolver,
               std::vector<std::vector<Complex>> basis,
               const std::vector<std::vector<Complex>>& portVectors);

    // 在指定频率评估 S 参数。
    SParameterPoint evaluate(double frequencyHz) const;

    // 在指定频率求解 reduced 坐标。
    std::vector<Complex> solveReduced(double frequencyHz) const;

    // 将 reduced 解重构为全阶边自由度。
    std::vector<Complex> reconstructField(double frequencyHz) const;

    // 返回 ROM 维度。
    int dimension() const { return romDim_; }

    // 返回原始全阶自由度维度。
    int fullDimension() const { return fullDim_; }

    // 判断 ROM 是否已经构造完成。
    bool ready() const { return ready_; }

    // 返回 ROM 基的最大正交性误差。
    double basisOrthogonalityError() const;

    // 返回最近一次 reduced solve 是否成功。
    bool lastSolveSucceeded() const { return lastSolveSucceeded_; }

private:
    const ProjectDefinition* project_ = nullptr;
    const FEMAssembler::AffineSystem* affine_ = nullptr;
    const PortModeSolver* portModeSolver_ = nullptr;

    bool ready_ = false;
    mutable bool lastSolveSucceeded_ = true;
    int romDim_ = 0;
    int fullDim_ = 0;
    int numProjectPorts_ = 0;
    int numVirtualPorts_ = 0;

    std::vector<std::vector<Complex>> basis_;
    std::vector<Complex> Ktilde_;
    std::vector<Complex> Mtilde_;
    std::vector<std::vector<Complex>> portModeReduced_;
};

}  // namespace fem::fastsweep
