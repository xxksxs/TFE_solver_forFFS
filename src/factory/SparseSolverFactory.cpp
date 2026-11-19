#include "bpfem/factory/SparseSolverFactory.hpp"

#include "bpfem/linalg/BiCGStabBackend.hpp"
#include "bpfem/linalg/GmresBackend.hpp"
#include "bpfem/linalg/MumpsBackend.hpp"
#include "bpfem/linalg/PardisoBackend.hpp"
#include "bpfem/linalg/PreconIlu0.hpp"
#include "bpfem/linalg/PreconJacobi.hpp"
#include "bpfem/linalg/PreconPardisoSingle.hpp"
#include "bpfem/linalg/PreconSchwarz.hpp"

#include <iostream>
#include <stdexcept>

namespace fem::factory {

namespace {

linalg::MumpsMode mumpsMode(const Options& options) {
    if (options.mumpsMode == "blr") return linalg::MumpsMode::BlockLowRank;
    if (options.mumpsMode == "ooc") return linalg::MumpsMode::OutOfCore;
    return linalg::MumpsMode::InCore;
}

std::shared_ptr<linalg::IPreconditioner> makePreconditioner(const Options& options) {
    switch (options.preconditioner) {
        case PreconditionerKind::Jacobi:
            return std::make_shared<linalg::PreconJacobi>();
        case PreconditionerKind::Ilu0:
            return std::make_shared<linalg::PreconIlu0>();
        case PreconditionerKind::PardisoSingle:
            return std::make_shared<linalg::PreconPardisoSingle>();
        case PreconditionerKind::Schwarz: {
            linalg::PreconSchwarz::LocalSolverFactory createLocalSolver;
            if (options.schwarzLocalMumps) {
#ifdef BPFEM_USE_MUMPS
                const auto mode = mumpsMode(options);
                if (mode == linalg::MumpsMode::OutOfCore) {
                    throw std::invalid_argument("Local MUMPS supports incore or blr; ooc is for global solves");
                }
                const auto tolerance = options.mumpsBlrTolerance;
                const auto path = options.outDir;
                createLocalSolver = [mode, tolerance, path] {
                    return std::make_unique<linalg::MumpsBackend>(mode, tolerance, path);
                };
#else
                throw std::invalid_argument("Schwarz with local MUMPS requires -DBPFEM_USE_MUMPS=ON");
#endif
            } else {
#ifdef BPFEM_USE_MKL
                createLocalSolver = [] { return std::make_unique<linalg::PardisoBackend>(); };
#else
                throw std::invalid_argument("Schwarz with local PARDISO requires Intel oneMKL");
#endif
            }
            return std::make_shared<linalg::PreconSchwarz>(
                options.schwarzSubdomains, options.schwarzOverlap, std::move(createLocalSolver));
        }
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

void warnDirectOptions(const Options& options, const char* backendName) {
    warnPreconIgnored(options.preconditioner, backendName);
    if (options.iterativeMklSpmv) {
        std::cerr << "[warn] --iterative-spmv mkl applies to GMRES/BiCGSTAB; "
                  << backendName << " uses its own sparse kernels\n";
    }
}

}  // namespace

LinearSolverKind resolveLinearSolverKind(LinearSolverKind requested) {
    if (requested != LinearSolverKind::Auto) return requested;
#ifdef BPFEM_USE_MKL
    return LinearSolverKind::Pardiso;
#else
    return LinearSolverKind::BiCGStab;
#endif
}

std::unique_ptr<linalg::ISparseSolver> makeSparseSolver(const Options& options) {
    const auto selected = resolveLinearSolverKind(options.linearSolver);
    if (options.preconditioner == PreconditionerKind::PardisoSingle
        && selected != LinearSolverKind::Gmres) {
        throw std::invalid_argument("--precon pardiso-single requires --linear-solver gmres");
    }
#ifndef BPFEM_USE_MKL
    if (options.iterativeMklSpmv) {
        throw std::invalid_argument("--iterative-spmv mkl requires a build with Intel oneMKL");
    }
#endif
    if (options.pardisoOutOfCore) {
#ifndef BPFEM_USE_MKL
        throw std::invalid_argument("--pardiso-mode ooc requires a build with Intel oneMKL");
#else
        if (selected != LinearSolverKind::Pardiso) {
            throw std::invalid_argument("--pardiso-mode ooc requires PARDISO (--linear-solver pardiso or auto)");
        }
#endif
    }
    auto buildBiCGStab = [&]() -> std::unique_ptr<linalg::ISparseSolver> {
        auto backend = std::make_unique<linalg::BiCGStabBackend>();
        backend->setMklSpmv(options.iterativeMklSpmv);
        if (auto precon = makePreconditioner(options); precon) {
            backend->setPreconditioner(std::move(precon));
        }
        return backend;
    };

    auto buildGmres = [&]() -> std::unique_ptr<linalg::ISparseSolver> {
        auto backend = std::make_unique<linalg::GmresBackend>();
        backend->setMklSpmv(options.iterativeMklSpmv);
        backend->setRestartLength(options.gmresRestart);
        if (auto precon = makePreconditioner(options); precon) {
            backend->setPreconditioner(std::move(precon));
        }
        return backend;
    };

    switch (selected) {
        case LinearSolverKind::Pardiso:
            warnDirectOptions(options, "PARDISO");
            return std::make_unique<linalg::PardisoBackend>(options.pardisoOutOfCore);
        case LinearSolverKind::BiCGStab:
            return buildBiCGStab();
        case LinearSolverKind::Gmres:
            return buildGmres();
        case LinearSolverKind::Mumps:
            warnDirectOptions(options, "MUMPS");
            return std::make_unique<linalg::MumpsBackend>(
                mumpsMode(options), options.mumpsBlrTolerance, options.outDir);
        case LinearSolverKind::Auto:
            break;
    }
    throw std::logic_error("Unresolved linear solver selection");
}

}  // namespace fem::factory
