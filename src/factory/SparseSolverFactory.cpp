#include "bpfem/factory/SparseSolverFactory.hpp"

#include "bpfem/linalg/BiCGStabBackend.hpp"
#include "bpfem/linalg/GmresBackend.hpp"
#include "bpfem/linalg/PardisoBackend.hpp"
#include "bpfem/linalg/PreconIlu0.hpp"
#include "bpfem/linalg/PreconJacobi.hpp"

#include <iostream>
#include <stdexcept>

namespace fem::factory {

namespace {

std::shared_ptr<linalg::IPreconditioner> makePreconditioner(PreconditionerKind kind) {
    switch (kind) {
        case PreconditionerKind::Jacobi:
            return std::make_shared<linalg::PreconJacobi>();
        case PreconditionerKind::Ilu0:
            return std::make_shared<linalg::PreconIlu0>();
        case PreconditionerKind::None:
        default:
            return nullptr;
    }
}

void warnPreconIgnored(PreconditionerKind kind, const char* backendName) {
    if (kind == PreconditionerKind::None) {
        return;
    }
    std::cerr << "[warn] --precon is ignored when --linear-solver selects '"
              << backendName << "' (preconditioning is iterative-solver only)\n";
}

}  // namespace

std::unique_ptr<linalg::ISparseSolver> makeSparseSolver(const Options& options) {
    auto buildBiCGStab = [&]() -> std::unique_ptr<linalg::ISparseSolver> {
        auto backend = std::make_unique<linalg::BiCGStabBackend>();
        if (auto precon = makePreconditioner(options.preconditioner); precon) {
            backend->setPreconditioner(std::move(precon));
        }
        return backend;
    };

    auto buildGmres = [&]() -> std::unique_ptr<linalg::ISparseSolver> {
        auto backend = std::make_unique<linalg::GmresBackend>();
        backend->setRestartLength(options.gmresRestart);
        if (auto precon = makePreconditioner(options.preconditioner); precon) {
            backend->setPreconditioner(std::move(precon));
        }
        return backend;
    };

    switch (options.linearSolver) {
        case LinearSolverKind::Direct:
            warnPreconIgnored(options.preconditioner, "PARDISO");
            return std::make_unique<linalg::PardisoBackend>();
        case LinearSolverKind::BiCGStab:
            return buildBiCGStab();
        case LinearSolverKind::Gmres:
            return buildGmres();
        case LinearSolverKind::Auto:
        default:
#ifdef BPFEM_USE_MKL
            warnPreconIgnored(options.preconditioner, "PARDISO");
            return std::make_unique<linalg::PardisoBackend>();
#else
            return buildBiCGStab();
#endif
    }
}

}  // namespace fem::factory
