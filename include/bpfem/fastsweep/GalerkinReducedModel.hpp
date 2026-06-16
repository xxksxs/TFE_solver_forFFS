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

    void build(const ProjectDefinition& project,
               const FEMAssembler::AffineSystem& affine,
               const PortModeSolver& portModeSolver,
               std::vector<std::vector<Complex>> basis,
               const std::vector<std::vector<Complex>>& portVectors);

    SParameterPoint evaluate(double frequencyHz) const;
    std::vector<Complex> solveReduced(double frequencyHz) const;
    std::vector<Complex> reconstructField(double frequencyHz) const;

    int dimension() const { return romDim_; }
    int fullDimension() const { return fullDim_; }
    bool ready() const { return ready_; }
    double basisOrthogonalityError() const;
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
