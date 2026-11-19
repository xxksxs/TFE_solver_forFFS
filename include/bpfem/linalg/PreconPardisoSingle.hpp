#pragma once

#include "bpfem/linalg/IPreconditioner.hpp"
#include "bpfem/linalg/MklPardisoSolver.hpp"

namespace fem::linalg {

// A single-precision factor used by flexible GMRES as a right preconditioner.
class PreconPardisoSingle final : public IPreconditioner {
public:
    PreconPardisoSingle();
    void setup(const SparseMatrix& matrix) override;
    void clear() noexcept override;
    void apply(const std::vector<std::complex<double>>& input,
               std::vector<std::complex<double>>& output) const override;
    const char* name() const override { return "PARDISO single"; }
    bool requiresFlexibleGmres() const override { return true; }

private:
    mutable MklPardisoSolver solver_{false, MklPardisoSolver::Precision::Single};
};

}  // namespace fem::linalg
