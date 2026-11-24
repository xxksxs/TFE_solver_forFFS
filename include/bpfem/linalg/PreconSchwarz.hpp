#pragma once

#include "bpfem/linalg/IFactorizedSparseSolver.hpp"
#include "bpfem/linalg/IPreconditioner.hpp"

#include <functional>
#include <memory>

namespace fem::linalg {

class PreconSchwarz final : public IPreconditioner {
public:
    using LocalSolverFactory = std::function<std::unique_ptr<IFactorizedSparseSolver>()>;

    PreconSchwarz(int subdomains, int overlap, LocalSolverFactory createLocalSolver);

    void setup(const SparseMatrix& matrix) override;
    void clear() noexcept override;
    void apply(const std::vector<std::complex<double>>& x,
               std::vector<std::complex<double>>& y) const override;
    const char* name() const override { return "Schwarz"; }

private:
    struct Subdomain {
        std::vector<int> indices;
        std::size_t coreSize = 0;
        std::unique_ptr<IFactorizedSparseSolver> solver;
    };

    int subdomainCount_;
    int overlap_;
    LocalSolverFactory createLocalSolver_;
    std::size_t dimension_ = 0;
    std::vector<Subdomain> subdomains_;
};

}  // namespace fem::linalg
