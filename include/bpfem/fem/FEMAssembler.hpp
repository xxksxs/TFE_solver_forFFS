#pragma once

#include "bpfem/bc/IBoundaryCondition.hpp"
#include "bpfem/core/Types.hpp"
#include "bpfem/fem/EdgeTopology.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#include <complex>
#include <map>
#include <memory>
#include <unordered_set>
#include <vector>

namespace fem {

class FEMAssembler {
public:
    FEMAssembler(const Mesh& mesh, const ProjectDefinition& project, const EdgeTopology& topology, const PortModeSolver& portModeSolver);

    // Register the additive boundary conditions (wave ports, ABCs, SIBCs, ...)
    // that contribute to the system matrix and RHS. PEC zero-Dirichlet is
    // *not* a boundary condition; it is collected automatically from
    // surface faces that are not waveports. Order of registration matches
    // the order BCs are applied (apply / declareSparsity / affineContributions).
    //
    // Default: no BCs registered. Application.cpp must register at least
    // a WavePortBC for waveport projects.
    void setBoundaryConditions(std::vector<std::shared_ptr<bc::IBoundaryCondition>> bcs);

    SparseMatrix assemble(double frequencyHz, std::vector<std::complex<double>>& rhs) const;

    // Frequency-independent affine decomposition of A(omega) for fast-sweep
    // model-order reduction (ALPS / Krylov MOR). For lossless problems
    //     A(omega) = K - k0(omega)^2 * M + j * sum_p beta_p(omega) * m_p m_p^T
    // where K and M are real-valued and pre-built from the element cache.
    // Conductive losses are not yet supported by this affine path; if any
    // material has a non-zero conductivity, AffineSystem::lossless is set to
    // false and the caller should fall back to the per-frequency direct path.
    struct AffineSystem {
        SparseMatrix K;                                                          // sum_e (1/mu_r,e) * K_e
        SparseMatrix M;                                                          // sum_e (eps_r,e) * M_e
        std::vector<std::vector<std::pair<int, double>>> portCoupling;           // m_p as sparse vectors (one per virtual port = (face, mode))
        std::vector<double> portCutoffSquared;                                   // k_{c,p}^2
        std::vector<int> portFaceIds;                                            // physical face id this virtual port belongs to
        std::vector<int> projectPortIndex;                                       // index into ProjectDefinition::ports (drives excitation phase / magnitude)
        std::vector<bool> isExcitationMode;                                      // true only for the excited mode of an excited port
        std::unordered_set<int> constrainedDofs;                                 // PEC zero-Dirichlet mask
        bool lossless = true;                                                    // false if any sigma > 0
    };

    AffineSystem buildAffineSystem() const;

private:
    // Per-element frequency-independent reference matrices. The volume integrand
    //     curlCurl(a,b) = (1/mu_r) * (curl N_a) . (curl N_b)
    //     mass    (a,b) =                   N_a   .       N_b
    // factors as
    //     K_e[a,b] = integral_T (curl N_a) . (curl N_b) dV       (geometry + basis only)
    //     M_e[a,b] = integral_T       N_a  .       N_b  dV       (geometry + basis only)
    // The frequency loop then assembles
    //     A_local(a,b)(omega) = invMu * K_e[a,b] - k0^2 * epsC(omega) * M_e[a,b]
    // which is O(localDof^2) instead of O(localDof^2 * Nq). Cached lazily on the
    // first assemble() call; sized once for the lifetime of the assembler since
    // mesh/topology are fixed.
    struct ElementMatrices {
        bool valid = false;     // false for degenerate or unmapped tets
        int dof = 0;            // localDof, mirrored from EdgeTopology
        double invMu = 1.0;     // 1 / relative permeability
        double epsR = 1.0;      // relative permittivity
        double conductivity = 0.0;  // S/m, multiplied by -j/(omega eps0) at solve time
        // Upper-triangle packed (size dof*(dof+1)/2) using upperTriIndex(a,b,dof).
        std::vector<double> curlCurlGeo;
        std::vector<double> massGeo;
    };

    const Mesh& mesh_;
    const ProjectDefinition& project_;
    const EdgeTopology& topology_;
    const PortModeSolver& portModeSolver_;
    std::vector<std::shared_ptr<bc::IBoundaryCondition>> bcs_;

    Material materialForBody(int bodyId) const;
    void applyBoundaryConditions(SparseMatrix& matrix, std::vector<std::complex<double>>& rhs, double frequencyHz, const std::unordered_set<int>& constrainedDofs) const;
    std::unordered_set<int> collectPecConstrainedDofs() const;
    std::shared_ptr<const SparsePattern> sparsityPattern(const std::unordered_set<int>& constrainedDofs) const;
    std::map<int, std::unordered_set<int>> collectSurfaceEdgesByFace() const;
    void ensureElementCache() const;
    const std::unordered_set<int>& cachedConstrainedDofs() const;

    SparseMatrix assembleVolumePiece(const std::shared_ptr<const SparsePattern>& pattern,
                                     const std::unordered_set<int>& constrainedDofs,
                                     bool useCurlCurl) const;
    std::shared_ptr<const SparsePattern> volumePatternFor(const std::unordered_set<int>& constrainedDofs) const;

    mutable std::shared_ptr<const SparsePattern> sparsityPattern_;
    mutable std::unordered_set<int> sparsityPatternConstrainedDofs_;
    mutable std::vector<ElementMatrices> elementCache_;
    mutable bool elementCacheBuilt_ = false;
    mutable std::unordered_set<int> constrainedDofsCache_;
    mutable bool constrainedDofsCached_ = false;
};

}  // namespace fem
