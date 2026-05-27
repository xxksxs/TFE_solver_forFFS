#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fem/EdgeTopology.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#include <complex>
#include <cstddef>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fem::bc {

// Per-call data passed to an IBoundaryCondition. Cheap to construct: holds
// references only.
//
// `constrainedDofs` is the PEC zero-Dirichlet mask collected by FEMAssembler.
// BCs MUST skip any (row, col) that lies in this set, both during assembly
// (`apply`) and pattern declaration (`declareSparsity`), so that the final
// imposeZeroDirichlet step does not have to undo work.
struct AssemblyContext {
    double frequencyHz = 0.0;
    const Mesh& mesh;
    const ProjectDefinition& project;
    const EdgeTopology& topology;
    const PortModeSolver& portModeSolver;
    const std::unordered_set<int>& constrainedDofs;
};

// Builder used during the sparsity-pattern phase. Wraps the per-row column
// sets that FEMAssembler::sparsityPattern is filling in. Each BC adds the
// (row, col) pairs its `apply` will later write to. Constrained rows or
// columns are silently dropped, matching FEMAssembler's existing rule.
class SparsePatternBuilder {
public:
    SparsePatternBuilder(std::size_t n,
                         std::vector<std::unordered_set<int>>& columnsByRow,
                         const std::unordered_set<int>& constrainedDofs)
        : n_(n), columnsByRow_(columnsByRow), constrainedDofs_(constrainedDofs) {}

    // Declares a non-zero at (row, col). The pair is stored upper-triangular
    // (smaller index in the row), matching SparseMatrix's storage. Pairs with
    // either index out of range or in the constrained set are silently
    // dropped, so callers can pass raw triangle DOFs without filtering.
    void add(int row, int col) {
        if (row < 0 || col < 0) return;
        if (static_cast<std::size_t>(row) >= n_) return;
        if (static_cast<std::size_t>(col) >= n_) return;
        if (constrainedDofs_.count(row) != 0) return;
        if (constrainedDofs_.count(col) != 0) return;
        if (col < row) std::swap(row, col);
        columnsByRow_[static_cast<std::size_t>(row)].insert(col);
    }

    std::size_t dimension() const { return n_; }

private:
    std::size_t n_;
    std::vector<std::unordered_set<int>>& columnsByRow_;
    const std::unordered_set<int>& constrainedDofs_;
};

// Sparse vector contribution that an additive port-style boundary makes to the
// frequency-independent affine decomposition consumed by ALPS. A single
// physical port may emit several of these (one per analytic mode in TFE
// multi-mode). Only one entry per project port should carry
// `isExcitationMode = true`.
struct AffinePortContribution {
    std::vector<std::pair<int, double>> coupling;  // m_p as (dof, weight)
    double cutoffSquared = 0.0;                    // k_{c,p}^2
    int faceId = -1;
    int projectPortIndex = -1;
    bool isExcitationMode = false;
};

// Abstract boundary condition that contributes additively to the FEM system.
//
// Concrete kinds (status):
//   - WavePortBC          implemented; replaces FEMAssembler::applyWavePorts
//   - AbsorbingBC         placeholder skeleton (1st-order ABC)
//   - ImpedanceBC         placeholder skeleton (Z_s surface impedance)
//   - FiniteConductorBC   placeholder skeleton (Leontovich SIBC)
//
// PEC zero-Dirichlet is *not* an IBoundaryCondition: it is a projection
// (rows/cols removed from the system) rather than an additive contribution,
// and the constrained-edge set is topology-driven rather than strategy-driven.
// PEC handling stays in FEMAssembler.
class IBoundaryCondition {
public:
    virtual ~IBoundaryCondition() = default;

    // Add this boundary's contribution to `matrix` and `rhs` for the given
    // frequency. Implementations must respect `ctx.constrainedDofs`.
    virtual void apply(SparseMatrix& matrix,
                       std::vector<std::complex<double>>& rhs,
                       const AssemblyContext& ctx) const = 0;

    // Declare every (row, col) pair `apply` will later write to, before the
    // CSR pattern is finalized. Skipping this step is allowed but kills
    // PARDISO's symbolic-factorization cache reuse across frequency points.
    virtual void declareSparsity(SparsePatternBuilder& builder,
                                 const AssemblyContext& ctx) const = 0;

    // For ALPS / fast-sweep MOR. Boundaries that decompose into rank-1 port
    // operators with frequency-independent shape report their port pieces
    // here; ABC / SIBC etc. that don't fit the affine form return an empty
    // vector. Default = empty.
    virtual std::vector<AffinePortContribution> affineContributions(
        const AssemblyContext& ctx) const {
        (void)ctx;
        return {};
    }

    // Human-readable name for logging.
    virtual const char* name() const = 0;
};

}  // namespace fem::bc

