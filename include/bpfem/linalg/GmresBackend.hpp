#pragma once

#include "bpfem/linalg/GmresSolver.hpp"
#include "bpfem/linalg/IPreconditioner.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"

#include <memory>

namespace fem::linalg {

// ISparseSolver wrapper around fem::GmresSolver. Composition keeps
// GmresSolver's signature stable.
class GmresBackend : public ISparseSolver {
public:
    GmresBackend() = default;

    void setPreconditioner(std::shared_ptr<IPreconditioner> precon);
    void setRestartLength(int m) { solver_.setRestartLength(m); }

    SolveResult solve(const SparseMatrix& A,
                      const std::vector<std::complex<double>>& b,
                      const SolverConfig& cfg = {}) override;

    const char* name() const override { return "GMRES"; }

private:
    GmresSolver solver_;
};

}  // namespace fem::linalg

