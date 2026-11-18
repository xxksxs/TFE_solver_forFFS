#pragma once

#include "bpfem/linalg/SparseMatrix.hpp"

#include <complex>
#include <vector>

namespace fem::linalg {

// Right-preconditioner contract: produce y = M^{-1} x.
//
// Iterative backends that support preconditioning call setup(A) once before
// the iteration loop and apply(x, y) inside the loop. Implementations may
// cache decomposition state from setup() across many apply() calls; the
// matrix passed to setup() must outlive subsequent apply() calls.
//
// Allowed alias: in-place apply (&x == &y) is supported by all in-tree
// implementations.
class IPreconditioner {
public:
    virtual ~IPreconditioner() = default;

    // Extract the data needed to apply M^{-1} from `A` (e.g. take the
    // diagonal for Jacobi; build an ILU factorization, etc). Called once per
    // solve. The reference must remain valid for subsequent apply() calls.
    virtual void setup(const SparseMatrix& A) = 0;
    virtual void clear() noexcept {}

    // GMRES must retain each M^{-1}v when apply() does not preserve a
    // linear combination to the arithmetic precision of the outer solver.
    virtual bool requiresFlexibleGmres() const { return false; }

    // y = M^{-1} x. y is resized to x.size() if needed.
    virtual void apply(const std::vector<std::complex<double>>& x,
                       std::vector<std::complex<double>>& y) const = 0;

    // Human-readable name for logging.
    virtual const char* name() const = 0;
};

}  // namespace fem::linalg

