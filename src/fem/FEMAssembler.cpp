#include "bpfem/fem/FEMAssembler.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/core/Math.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace fem {

namespace {

struct BasisValue {
    Vec3 value;
    Vec3 curl;
};

struct QuadraturePoint {
    std::array<double, 4> lambda;
    double weight = 0.0;
};

std::array<QuadraturePoint, 11> tetraQuadraturePoints() {
    constexpr double a = 0.7857142857142857;
    constexpr double b = 0.0714285714285714;
    constexpr double c = 0.3994035761667992;
    constexpr double d = 0.1005964238332008;
    return {{{{{0.25, 0.25, 0.25, 0.25}}, -0.0789333333333333},
             {{{a, b, b, b}}, 0.0457333333333333},
             {{{b, a, b, b}}, 0.0457333333333333},
             {{{b, b, a, b}}, 0.0457333333333333},
             {{{b, b, b, a}}, 0.0457333333333333},
             {{{c, c, d, d}}, 0.1493333333333333},
             {{{c, d, c, d}}, 0.1493333333333333},
             {{{c, d, d, c}}, 0.1493333333333333},
             {{{d, c, c, d}}, 0.1493333333333333},
             {{{d, c, d, c}}, 0.1493333333333333},
             {{{d, d, c, c}}, 0.1493333333333333}}};
}

BasisValue evaluateBasis(const LocalDofRef& ref, const std::array<double, 4>& lambda, const std::array<Vec3, 4>& grad) {
    const int a = ref.localNodes[0];
    const int b = ref.localNodes[1];
    const int c = ref.localNodes[2];
    if (a < 0 || b < 0) {
        return {};
    }

    const Vec3 zeroValue = lambda[static_cast<std::size_t>(a)] * grad[static_cast<std::size_t>(b)]
                         - lambda[static_cast<std::size_t>(b)] * grad[static_cast<std::size_t>(a)];
    const Vec3 zeroCurl = 2.0 * cross(grad[static_cast<std::size_t>(a)], grad[static_cast<std::size_t>(b)]);

    if (ref.kind == LocalDofKind::EdgeZero) {
        return {ref.sign * zeroValue, ref.sign * zeroCurl};
    }
    if (ref.kind == LocalDofKind::EdgeFirst) {
        const double factor = lambda[static_cast<std::size_t>(a)] - lambda[static_cast<std::size_t>(b)];
        const Vec3 factorGrad = grad[static_cast<std::size_t>(a)] - grad[static_cast<std::size_t>(b)];
        return {factor * zeroValue, cross(factorGrad, zeroValue) + factor * zeroCurl};
    }
    if (c < 0) {
        return {};
    }

    const Vec3 bcValue = lambda[static_cast<std::size_t>(b)] * grad[static_cast<std::size_t>(c)]
                       - lambda[static_cast<std::size_t>(c)] * grad[static_cast<std::size_t>(b)];
    const Vec3 bcCurl = 2.0 * cross(grad[static_cast<std::size_t>(b)], grad[static_cast<std::size_t>(c)]);
    if (ref.kind == LocalDofKind::FaceFirst0) {
        const double factor = lambda[static_cast<std::size_t>(c)];
        const Vec3 factorGrad = grad[static_cast<std::size_t>(c)];
        return {factor * zeroValue, cross(factorGrad, zeroValue) + factor * zeroCurl};
    }

    // FaceFirst1: lambda_a * N_bc, linearly independent from FaceFirst0
    // (lambda_c * N_ab). The previous form lambda_a*N_bc - lambda_b*N_ac
    // collapses to -lambda_c * N_ab = -FaceFirst0, leaving the local face
    // bubble pair rank-deficient. PortModeSolver carried the same mistake.
    const double factor = lambda[static_cast<std::size_t>(a)];
    const Vec3 factorGrad = grad[static_cast<std::size_t>(a)];
    return {factor * bcValue, cross(factorGrad, bcValue) + factor * bcCurl};
}

// Upper-triangle packed index for an n x n symmetric matrix:
//   (a, b) with a <= b -> a * (2*n - a - 1) / 2 + b
// Storage size = n*(n+1)/2.
inline std::size_t upperTriIndex(int a, int b, int n) {
    if (a > b) {
        std::swap(a, b);
    }
    return static_cast<std::size_t>(a) * static_cast<std::size_t>(2 * n - a - 1) / 2u
         + static_cast<std::size_t>(b);
}

}  // namespace

FEMAssembler::FEMAssembler(const Mesh& mesh, const ProjectDefinition& project, const EdgeTopology& topology, const PortModeSolver& portModeSolver)
    : mesh_(mesh), project_(project), topology_(topology), portModeSolver_(portModeSolver) {}

void FEMAssembler::setBoundaryConditions(std::vector<std::shared_ptr<bc::IBoundaryCondition>> bcs) {
    bcs_ = std::move(bcs);
    // Sparsity pattern depends on the BC set, so invalidate the cached one.
    sparsityPattern_.reset();
    sparsityPatternConstrainedDofs_.clear();
}

void FEMAssembler::ensureElementCache() const {
    if (elementCacheBuilt_) {
        return;
    }
    const auto& elementDofs = topology_.elementDofs();
    const auto quadraturePoints = tetraQuadraturePoints();
    elementCache_.assign(mesh_.tetrahedra.size(), ElementMatrices{});

    for (std::size_t elementIndex = 0; elementIndex < mesh_.tetrahedra.size(); ++elementIndex) {
        const auto& tet = mesh_.tetrahedra[elementIndex];
        std::array<Vec3, 4> p{};
        bool valid = true;
        for (int i = 0; i < 4; ++i) {
            const auto pit = mesh_.pointsById.find(tet.vertexIds[static_cast<std::size_t>(i)]);
            if (pit == mesh_.pointsById.end()) {
                valid = false;
                break;
            }
            p[static_cast<std::size_t>(i)] = pit->second;
        }
        if (!valid) {
            continue;
        }
        const double volume = tetraVolume(p[0], p[1], p[2], p[3]);
        if (volume <= 1.0e-24) {
            continue;
        }

        const auto grad = tetraGradients(p);
        const auto& refs = elementDofs[elementIndex];
        const int dof = static_cast<int>(refs.size());
        const std::size_t packedSize = static_cast<std::size_t>(dof) * static_cast<std::size_t>(dof + 1) / 2u;
        ElementMatrices entry;
        entry.valid = true;
        entry.dof = dof;
        entry.curlCurlGeo.assign(packedSize, 0.0);
        entry.massGeo.assign(packedSize, 0.0);
        const Material material = materialForBody(tet.bodyId);
        entry.invMu = 1.0 / std::max(material.relativePermeability, 1.0e-30);
        entry.epsR = material.relativePermittivity;
        entry.conductivity = material.conductivity;

        for (int a = 0; a < dof; ++a) {
            for (int b = a; b < dof; ++b) {
                double curlCurl = 0.0;
                double mass = 0.0;
                for (const auto& qp : quadraturePoints) {
                    const BasisValue rowBasis = evaluateBasis(refs[static_cast<std::size_t>(a)], qp.lambda, grad);
                    const BasisValue colBasis = evaluateBasis(refs[static_cast<std::size_t>(b)], qp.lambda, grad);
                    curlCurl += qp.weight * volume * dot(rowBasis.curl, colBasis.curl);
                    mass     += qp.weight * volume * dot(rowBasis.value, colBasis.value);
                }
                const std::size_t idx = upperTriIndex(a, b, dof);
                entry.curlCurlGeo[idx] = curlCurl;
                entry.massGeo[idx] = mass;
            }
        }
        elementCache_[elementIndex] = std::move(entry);
    }
    elementCacheBuilt_ = true;
}

SparseMatrix FEMAssembler::assemble(double frequencyHz, std::vector<std::complex<double>>& rhs) const {
    ensureElementCache();
    const auto& constrainedDofs = cachedConstrainedDofs();
    SparseMatrix matrix(sparsityPattern(constrainedDofs));
    rhs.assign(topology_.edgeCount(), 0.0);
    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double k0Sq = k0 * k0;
    const double omegaEps0 = 2.0 * pi * frequencyHz * epsilon0;
    const auto& elementDofs = topology_.elementDofs();
    const auto isConstrained = [&constrainedDofs](int dof) {
        return constrainedDofs.count(dof) != 0;
    };

    for (std::size_t elementIndex = 0; elementIndex < mesh_.tetrahedra.size(); ++elementIndex) {
        const auto& entry = elementCache_[elementIndex];
        if (!entry.valid) {
            continue;
        }
        // Frequency-dependent material scalars; geometry-only K_e and M_e were
        // built once in ensureElementCache().
        const std::complex<double> epsC(entry.epsR, -entry.conductivity / std::max(omegaEps0, 1.0e-300));
        const double curlScale = entry.invMu;
        const std::complex<double> massScale = -k0Sq * epsC;
        const auto& refs = elementDofs[elementIndex];
        const int dof = entry.dof;

        for (int a = 0; a < dof; ++a) {
            const auto& row = refs[static_cast<std::size_t>(a)];
            if (row.globalIndex < 0 || isConstrained(row.globalIndex)) {
                continue;
            }
            for (int b = a; b < dof; ++b) {
                const auto& col = refs[static_cast<std::size_t>(b)];
                if (col.globalIndex < 0 || isConstrained(col.globalIndex)) {
                    continue;
                }
                const std::size_t idx = upperTriIndex(a, b, dof);
                const std::complex<double> contribution =
                    curlScale * entry.curlCurlGeo[idx] + massScale * entry.massGeo[idx];
                matrix.add(row.globalIndex, col.globalIndex, contribution);
            }
        }
    }

    applyBoundaryConditions(matrix, rhs, frequencyHz, constrainedDofs);
    matrix.imposeZeroDirichlet(constrainedDofs, rhs);
    return matrix;
}

FEMAssembler::AffineSystem FEMAssembler::buildAffineSystem() const {
    ensureElementCache();
    AffineSystem out;
    out.constrainedDofs = cachedConstrainedDofs();
    const auto pattern = volumePatternFor(out.constrainedDofs);
    out.K = assembleVolumePiece(pattern, out.constrainedDofs, true);
    out.M = assembleVolumePiece(pattern, out.constrainedDofs, false);

    // Detect lossy materials: if any tet has sigma > 0, the affine
    // decomposition with a single real M no longer captures the full
    // frequency dependence (epsC = epsR - j sigma/(omega eps0)).
    out.lossless = true;
    for (const auto& entry : elementCache_) {
        if (entry.valid && entry.conductivity > 0.0) {
            out.lossless = false;
            break;
        }
    }

    // Port coupling vectors and cutoff k_{c,p}^2 are now gathered through
    // the registered IBoundaryCondition instances. Each BC may emit zero or
    // more AffinePortContribution entries, in BC registration order. For TFE
    // multi-mode ports each (face, analytic mode) yields one entry; only the
    // excitation mode of an excited project port has isExcitationMode = true.
    bc::AssemblyContext ctx{
        0.0,
        mesh_,
        project_,
        topology_,
        portModeSolver_,
        out.constrainedDofs
    };
    std::size_t reserveHint = project_.ports.size();
    out.portFaceIds.reserve(reserveHint);
    out.portCoupling.reserve(reserveHint);
    out.portCutoffSquared.reserve(reserveHint);
    out.projectPortIndex.reserve(reserveHint);
    out.isExcitationMode.reserve(reserveHint);
    for (const auto& cond : bcs_) {
        if (!cond) continue;
        for (auto& contrib : cond->affineContributions(ctx)) {
            out.portFaceIds.push_back(contrib.faceId);
            out.portCutoffSquared.push_back(contrib.cutoffSquared);
            out.portCoupling.push_back(std::move(contrib.coupling));
            out.projectPortIndex.push_back(contrib.projectPortIndex);
            out.isExcitationMode.push_back(contrib.isExcitationMode);
        }
    }

    return out;
}

std::shared_ptr<const SparsePattern> FEMAssembler::volumePatternFor(const std::unordered_set<int>& constrainedDofs) const {
    // Rebuild a sparsity pattern that contains only volume-element couplings,
    // without the rank-1 port outer products. ALPS K and M never touch the
    // port term, so this pattern is strictly smaller (or equal) to the one
    // used by the per-frequency direct path.
    const auto isConstrained = [&constrainedDofs](int dof) {
        return constrainedDofs.count(dof) != 0;
    };
    const auto n = topology_.edgeCount();
    std::vector<std::unordered_set<int>> columnsByRowSet(n);
    for (const auto& refs : topology_.elementDofs()) {
        for (std::size_t a = 0; a < refs.size(); ++a) {
            for (std::size_t b = a; b < refs.size(); ++b) {
                int row = refs[a].globalIndex;
                int col = refs[b].globalIndex;
                if (row < 0 || col < 0 || static_cast<std::size_t>(row) >= n
                    || static_cast<std::size_t>(col) >= n
                    || isConstrained(row) || isConstrained(col)) {
                    continue;
                }
                if (col < row) {
                    std::swap(row, col);
                }
                columnsByRowSet[static_cast<std::size_t>(row)].insert(col);
            }
        }
    }
    for (int dof : constrainedDofs) {
        if (dof >= 0 && static_cast<std::size_t>(dof) < n) {
            columnsByRowSet[static_cast<std::size_t>(dof)].insert(dof);
        }
    }
    std::vector<std::vector<int>> columnsByRow(n);
    for (std::size_t row = 0; row < n; ++row) {
        columnsByRow[row].assign(columnsByRowSet[row].begin(), columnsByRowSet[row].end());
        std::sort(columnsByRow[row].begin(), columnsByRow[row].end());
    }
    return SparseMatrix::createPattern(n, columnsByRow);
}

SparseMatrix FEMAssembler::assembleVolumePiece(const std::shared_ptr<const SparsePattern>& pattern,
                                               const std::unordered_set<int>& constrainedDofs,
                                               bool useCurlCurl) const {
    SparseMatrix matrix(pattern);
    const auto& elementDofs = topology_.elementDofs();
    const auto isConstrained = [&constrainedDofs](int dof) {
        return constrainedDofs.count(dof) != 0;
    };
    for (std::size_t elementIndex = 0; elementIndex < mesh_.tetrahedra.size(); ++elementIndex) {
        const auto& entry = elementCache_[elementIndex];
        if (!entry.valid) {
            continue;
        }
        const double scale = useCurlCurl ? entry.invMu : entry.epsR;
        const auto& src = useCurlCurl ? entry.curlCurlGeo : entry.massGeo;
        const auto& refs = elementDofs[elementIndex];
        const int dof = entry.dof;
        for (int a = 0; a < dof; ++a) {
            const auto& row = refs[static_cast<std::size_t>(a)];
            if (row.globalIndex < 0 || isConstrained(row.globalIndex)) {
                continue;
            }
            for (int b = a; b < dof; ++b) {
                const auto& col = refs[static_cast<std::size_t>(b)];
                if (col.globalIndex < 0 || isConstrained(col.globalIndex)) {
                    continue;
                }
                const std::size_t idx = static_cast<std::size_t>(a) * static_cast<std::size_t>(2 * dof - a - 1) / 2u
                                       + static_cast<std::size_t>(b);
                matrix.add(row.globalIndex, col.globalIndex, std::complex<double>(scale * src[idx], 0.0));
            }
        }
    }
    return matrix;
}

Material FEMAssembler::materialForBody(int bodyId) const {
    auto bit = mesh_.bodies.find(bodyId);
    if (bit != mesh_.bodies.end()) {
        auto mit = project_.materials.find(bit->second.name);
        if (mit != project_.materials.end()) {
            return mit->second;
        }
    }
    auto mit = project_.materials.find(project_.backgroundMaterial);
    if (mit != project_.materials.end()) {
        return mit->second;
    }
    return Material{"vacuum", 1.0, 1.0, 0.0};
}

// Walk the registered IBoundaryCondition list and let each BC add its own
// contribution to the matrix and RHS. With the empty default registration
// nothing happens; FEMAssembler is therefore neutral with respect to
// boundary types and any new BC (waveports, ABC, SIBC, ...) is added in
// Application::runApplication via setBoundaryConditions().
//
// The PEC zero-Dirichlet projection runs *after* this method returns
// (assemble() calls SparseMatrix::imposeZeroDirichlet); BCs must skip
// (row, col) pairs that are in `constrainedDofs` themselves to avoid
// fighting with that projection.
void FEMAssembler::applyBoundaryConditions(SparseMatrix& matrix, std::vector<std::complex<double>>& rhs, double frequencyHz, const std::unordered_set<int>& constrainedDofs) const {
    if (bcs_.empty()) {
        return;
    }
    bc::AssemblyContext ctx{
        frequencyHz,
        mesh_,
        project_,
        topology_,
        portModeSolver_,
        constrainedDofs
    };
    for (const auto& cond : bcs_) {
        if (cond) {
            cond->apply(matrix, rhs, ctx);
        }
    }
}

std::unordered_set<int> FEMAssembler::collectPecConstrainedDofs() const {
    std::unordered_set<int> portFaces;
    std::unordered_set<int> constrainedEdges;
    for (const auto& port : project_.ports) {
        portFaces.insert(port.faceId);
    }
    for (const auto& [faceId, bodyEdges] : collectSurfaceEdgesByFace()) {
        if (portFaces.count(faceId) != 0) {
            continue;
        }
        constrainedEdges.insert(bodyEdges.begin(), bodyEdges.end());
    }
    return constrainedEdges;
}

const std::unordered_set<int>& FEMAssembler::cachedConstrainedDofs() const {
    if (!constrainedDofsCached_) {
        constrainedDofsCache_ = collectPecConstrainedDofs();
        constrainedDofsCached_ = true;
    }
    return constrainedDofsCache_;
}

std::shared_ptr<const SparsePattern> FEMAssembler::sparsityPattern(const std::unordered_set<int>& constrainedDofs) const {
    if (sparsityPattern_ && sparsityPatternConstrainedDofs_ == constrainedDofs) {
        return sparsityPattern_;
    }

    const auto isConstrained = [&constrainedDofs](int dof) {
        return constrainedDofs.count(dof) != 0;
    };
    const auto n = topology_.edgeCount();
    std::vector<std::unordered_set<int>> columnsByRowSet(n);
    auto addPatternEntry = [&](int row, int col) {
        if (row < 0 || col < 0 || static_cast<std::size_t>(row) >= n || static_cast<std::size_t>(col) >= n) {
            return;
        }
        if (isConstrained(row) || isConstrained(col)) {
            return;
        }
        if (col < row) {
            std::swap(row, col);
        }
        columnsByRowSet[static_cast<std::size_t>(row)].insert(col);
    };

    for (const auto& refs : topology_.elementDofs()) {
        for (std::size_t a = 0; a < refs.size(); ++a) {
            for (std::size_t b = a; b < refs.size(); ++b) {
                addPatternEntry(refs[a].globalIndex, refs[b].globalIndex);
            }
        }
    }

    if (!bcs_.empty()) {
        bc::AssemblyContext ctx{
            0.0,
            mesh_,
            project_,
            topology_,
            portModeSolver_,
            constrainedDofs
        };
        bc::SparsePatternBuilder builder(n, columnsByRowSet, constrainedDofs);
        for (const auto& cond : bcs_) {
            if (cond) {
                cond->declareSparsity(builder, ctx);
            }
        }
    }

    for (int dof : constrainedDofs) {
        if (dof >= 0 && static_cast<std::size_t>(dof) < n) {
            columnsByRowSet[static_cast<std::size_t>(dof)].insert(dof);
        }
    }

    std::vector<std::vector<int>> columnsByRow(n);
    for (std::size_t row = 0; row < n; ++row) {
        columnsByRow[row].assign(columnsByRowSet[row].begin(), columnsByRowSet[row].end());
        std::sort(columnsByRow[row].begin(), columnsByRow[row].end());
    }

    sparsityPattern_ = SparseMatrix::createPattern(n, columnsByRow);
    sparsityPatternConstrainedDofs_ = constrainedDofs;
    return sparsityPattern_;
}

std::map<int, std::unordered_set<int>> FEMAssembler::collectSurfaceEdgesByFace() const {
    std::map<int, std::unordered_set<int>> byFace;
    for (const auto& port : project_.ports) {
        byFace[port.faceId] = topology_.edgesOnFace(port.faceId);
    }
    for (const auto& tri : mesh_.surfaceTriangles) {
        const auto fit = mesh_.facets.find(tri.facetId);
        if (fit == mesh_.facets.end()) {
            continue;
        }
        const int faceId = fit->second.faceId;
        const auto faceEdges = topology_.edgesOnFace(faceId);
        byFace[faceId].insert(faceEdges.begin(), faceEdges.end());
    }
    return byFace;
}

}  // namespace fem
