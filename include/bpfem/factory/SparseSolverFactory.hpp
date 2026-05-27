#pragma once

#include "bpfem/app/Application.hpp"   // Options, LinearSolverKind
#include "bpfem/linalg/ISparseSolver.hpp"

#include <memory>

namespace fem::factory {

// Builds the sparse solver backend that matches Options::linearSolver.
//
// Routing rules:
//   - Options::linearSolver == LinearSolverKind::Auto (default):
//       prefers PARDISO when this binary is built with BPFEM_USE_MKL,
//       otherwise falls back to BiCGSTAB. Matches the legacy compile-time
//       behavior so unchanged invocations keep their old solver choice.
//   - LinearSolverKind::Direct:
//       returns PARDISO; throws at solve() time if BPFEM_USE_MKL is off.
//   - LinearSolverKind::BiCGStab:
//       returns BiCGSTAB unconditionally (forces iterative even on MKL builds).
//
// The factory function is the *only* place in the code base that should
// contain BPFEM_USE_MKL preprocessor branches related to solver selection.
std::unique_ptr<linalg::ISparseSolver> makeSparseSolver(const Options& options);

}  // namespace fem::factory

