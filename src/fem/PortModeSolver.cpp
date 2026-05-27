#include "bpfem/fem/PortModeSolver.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/core/Math.hpp"

#ifdef BPFEM_USE_MKL
#include <mkl.h>
#include <mkl_lapacke.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace fem {
namespace {

using Matrix = std::vector<std::vector<double>>;

struct PortDofRef {
    int portIndex = -1;
    int globalIndex = -1;
    LocalDofKind kind = LocalDofKind::EdgeZero;
    std::array<int, 3> localNodes{{-1, -1, -1}};
    double sign = 1.0;
};

struct BasisValue {
    Vec3 value;
    double curl = 0.0;
};

struct TriangleQuadraturePoint {
    std::array<double, 3> lambda;
    double weight = 0.0;  // normalized so that sum_i weight_i = 1
};

// Wandzura 7-point rule on the reference triangle, exact for polynomials of
// total degree <= 5. Sufficient for 1st-order hierarchical H(curl) port mass
// matrices, whose integrand reaches total degree 4 in barycentric coordinates.
// The 3-point degree-2 rule used previously under-integrated the mass matrix,
// causing LAPACKE_dsygv to report a non-positive-definite mass operator
// (info > N) on order-1 runs.
std::array<TriangleQuadraturePoint, 7> triangleQuadraturePoints() {
    constexpr double a1 = 0.0597158717897698;
    constexpr double b1 = 0.4701420641051151;
    constexpr double a2 = 0.7974269853530873;
    constexpr double b2 = 0.1012865073234563;
    constexpr double w0 = 0.225;
    constexpr double w1 = 0.13239415278850618;
    constexpr double w2 = 0.12593918054482716;
    return {{{{{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}}, w0},
             {{{a1, b1, b1}}, w1},
             {{{b1, a1, b1}}, w1},
             {{{b1, b1, a1}}, w1},
             {{{a2, b2, b2}}, w2},
             {{{b2, a2, b2}}, w2},
             {{{b2, b2, a2}}, w2}}};
}

Matrix makeMatrix(std::size_t n) {
    return Matrix(n, std::vector<double>(n, 0.0));
}

std::vector<double> multiply(const Matrix& a, const std::vector<double>& x) {
    std::vector<double> y(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < x.size(); ++j) {
            y[i] += a[i][j] * x[j];
        }
    }
    return y;
}

Matrix cholesky(Matrix m) {
    const std::size_t n = m.size();
    for (std::size_t i = 0; i < n; ++i) {
        m[i][i] += 1.0e-18;
    }
    Matrix l = makeMatrix(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = m[i][j];
            for (std::size_t k = 0; k < j; ++k) {
                sum -= l[i][k] * l[j][k];
            }
            if (i == j) {
                l[i][j] = std::sqrt(std::max(sum, 1.0e-30));
            } else {
                l[i][j] = sum / std::max(l[j][j], 1.0e-30);
            }
        }
    }
    return l;
}

std::vector<double> solveLower(const Matrix& l, const std::vector<double>& b) {
    std::vector<double> x(b.size(), 0.0);
    for (std::size_t i = 0; i < b.size(); ++i) {
        double sum = b[i];
        for (std::size_t j = 0; j < i; ++j) {
            sum -= l[i][j] * x[j];
        }
        x[i] = sum / std::max(l[i][i], 1.0e-30);
    }
    return x;
}

std::vector<double> solveUpper(const Matrix& l, const std::vector<double>& b) {
    std::vector<double> x(b.size(), 0.0);
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(b.size()) - 1; i >= 0; --i) {
        double sum = b[static_cast<std::size_t>(i)];
        for (std::size_t j = static_cast<std::size_t>(i) + 1; j < b.size(); ++j) {
            sum -= l[j][static_cast<std::size_t>(i)] * x[j];
        }
        x[static_cast<std::size_t>(i)] = sum / std::max(l[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)], 1.0e-30);
    }
    return x;
}

Matrix generalizedToStandard(const Matrix& k, const Matrix& m) {
    const std::size_t n = k.size();
    const Matrix l = cholesky(m);
    Matrix a = makeMatrix(n);
    for (std::size_t j = 0; j < n; ++j) {
        std::vector<double> e(n, 0.0);
        e[j] = 1.0;
        const auto invLtColumn = solveUpper(l, e);
        const auto kInvLtColumn = multiply(k, invLtColumn);
        const auto column = solveLower(l, kInvLtColumn);
        for (std::size_t i = 0; i < n; ++i) {
            a[i][j] = column[i];
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double value = 0.5 * (a[i][j] + a[j][i]);
            a[i][j] = value;
            a[j][i] = value;
        }
    }
    return a;
}

std::pair<std::vector<double>, Matrix> jacobiEigen(Matrix a) {
    const std::size_t n = a.size();
    Matrix v = makeMatrix(n);
    for (std::size_t i = 0; i < n; ++i) {
        v[i][i] = 1.0;
    }
    const int maxIterations = static_cast<int>(std::max<std::size_t>(100, 50 * n * n));
    for (int iter = 0; iter < maxIterations; ++iter) {
        std::size_t p = 0;
        std::size_t q = 1;
        double maxOffDiagonal = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                const double value = std::abs(a[i][j]);
                if (value > maxOffDiagonal) {
                    maxOffDiagonal = value;
                    p = i;
                    q = j;
                }
            }
        }
        if (maxOffDiagonal < 1.0e-10) {
            break;
        }
        const double theta = 0.5 * std::atan2(2.0 * a[p][q], a[q][q] - a[p][p]);
        const double c = std::cos(theta);
        const double s = std::sin(theta);
        const double app = a[p][p];
        const double aqq = a[q][q];
        const double apq = a[p][q];
        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
        for (std::size_t r = 0; r < n; ++r) {
            if (r == p || r == q) {
                continue;
            }
            const double arp = a[r][p];
            const double arq = a[r][q];
            a[r][p] = c * arp - s * arq;
            a[p][r] = a[r][p];
            a[r][q] = s * arp + c * arq;
            a[q][r] = a[r][q];
        }
        for (std::size_t r = 0; r < n; ++r) {
            const double vrp = v[r][p];
            const double vrq = v[r][q];
            v[r][p] = c * vrp - s * vrq;
            v[r][q] = s * vrp + c * vrq;
        }
    }
    std::vector<double> eigenvalues(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        eigenvalues[i] = a[i][i];
    }
    return {eigenvalues, v};
}

std::pair<std::vector<double>, Matrix> generalizedEigen(const Matrix& stiffness, const Matrix& mass) {
#ifdef BPFEM_USE_MKL
    const std::size_t n = stiffness.size();
    std::vector<double> a(n * n, 0.0);
    std::vector<double> b(n * n, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            a[col * n + row] = stiffness[row][col];
            b[col * n + row] = mass[row][col];
        }
    }
    std::vector<double> eigenvalues(n, 0.0);
    const auto order = static_cast<lapack_int>(n);
    const lapack_int info = LAPACKE_dsygv(LAPACK_COL_MAJOR, 1, 'V', 'U', order, a.data(), order, b.data(), order, eigenvalues.data());
    if (info != 0) {
        throw std::runtime_error("Port H(curl) eigenproblem failed with LAPACKE_dsygv info " + std::to_string(info));
    }
    Matrix eigenvectors = makeMatrix(n);
    for (std::size_t col = 0; col < n; ++col) {
        for (std::size_t row = 0; row < n; ++row) {
            eigenvectors[row][col] = a[col * n + row];
        }
    }
    return {eigenvalues, eigenvectors};
#else
    const Matrix standard = generalizedToStandard(stiffness, mass);
    const auto [eigenvalues, eigenvectorsStandard] = jacobiEigen(standard);
    const Matrix lowerMass = cholesky(mass);
    Matrix eigenvectors = makeMatrix(stiffness.size());
    for (std::size_t col = 0; col < stiffness.size(); ++col) {
        std::vector<double> y(stiffness.size(), 0.0);
        for (std::size_t row = 0; row < stiffness.size(); ++row) {
            y[row] = eigenvectorsStandard[row][col];
        }
        const auto dofs = solveUpper(lowerMass, y);
        for (std::size_t row = 0; row < stiffness.size(); ++row) {
            eigenvectors[row][col] = dofs[row];
        }
    }
    return {eigenvalues, eigenvectors};
#endif
}

std::array<Vec3, 3> surfaceGradients(const std::array<Vec3, 3>& p, double& area, Vec3& normal) {
    const Vec3 e01 = p[1] - p[0];
    const Vec3 e02 = p[2] - p[0];
    normal = cross(e01, e02);
    const double normalLength = norm(normal);
    area = 0.5 * normalLength;
    if (normalLength <= 1.0e-30 || norm(e01) <= 1.0e-30) {
        return {};
    }
    normal = normal / normalLength;
    const Vec3 u = e01 / norm(e01);
    const Vec3 v = cross(normal, u);
    const std::array<double, 3> x{0.0, norm(e01), dot(e02, u)};
    const std::array<double, 3> y{0.0, 0.0, dot(e02, v)};
    const double twoArea = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]);
    if (std::abs(twoArea) <= 1.0e-30) {
        area = 0.0;
        return {};
    }
    std::array<Vec3, 3> grad{};
    const std::array<double, 3> gx{(y[1] - y[2]) / twoArea, (y[2] - y[0]) / twoArea, (y[0] - y[1]) / twoArea};
    const std::array<double, 3> gy{(x[2] - x[1]) / twoArea, (x[0] - x[2]) / twoArea, (x[1] - x[0]) / twoArea};
    for (std::size_t i = 0; i < 3; ++i) {
        grad[i] = gx[i] * u + gy[i] * v;
    }
    area = std::abs(twoArea) * 0.5;
    return grad;
}

Vec3 edgeBasisValue(int a, int b, const std::array<double, 3>& lambda, const std::array<Vec3, 3>& grad) {
    return lambda[static_cast<std::size_t>(a)] * grad[static_cast<std::size_t>(b)]
         - lambda[static_cast<std::size_t>(b)] * grad[static_cast<std::size_t>(a)];
}

Vec3 edgeFirstBasisValue(int a, int b, const std::array<double, 3>& lambda, const std::array<Vec3, 3>& grad) {
    return (lambda[static_cast<std::size_t>(a)] - lambda[static_cast<std::size_t>(b)]) * edgeBasisValue(a, b, lambda, grad);
}

Vec3 faceBasisValue(int index, const std::array<int, 3>& nodes, const std::array<double, 3>& lambda, const std::array<Vec3, 3>& grad) {
    const int a = nodes[0];
    const int b = nodes[1];
    const int c = nodes[2];
    if (index == 0) {
        return lambda[static_cast<std::size_t>(c)] * edgeBasisValue(a, b, lambda, grad);
    }
    return lambda[static_cast<std::size_t>(a)] * edgeBasisValue(b, c, lambda, grad)
         - lambda[static_cast<std::size_t>(b)] * edgeBasisValue(a, c, lambda, grad);
}

BasisValue evaluateBasis(const PortDofRef& ref, const std::array<double, 3>& lambda, const std::array<Vec3, 3>& grad, const Vec3& normal) {
    const int a = ref.localNodes[0];
    const int b = ref.localNodes[1];
    const int c = ref.localNodes[2];
    if (a < 0 || b < 0) {
        return {};
    }
    const Vec3 zeroValue = edgeBasisValue(a, b, lambda, grad);
    const double zeroCurl = 2.0 * dot(normal, cross(grad[static_cast<std::size_t>(a)], grad[static_cast<std::size_t>(b)]));
    if (ref.kind == LocalDofKind::EdgeZero) {
        return {ref.sign * zeroValue, ref.sign * zeroCurl};
    }
    if (ref.kind == LocalDofKind::EdgeFirst) {
        const double factor = lambda[static_cast<std::size_t>(a)] - lambda[static_cast<std::size_t>(b)];
        const Vec3 factorGrad = grad[static_cast<std::size_t>(a)] - grad[static_cast<std::size_t>(b)];
        return {factor * zeroValue, dot(normal, cross(factorGrad, zeroValue)) + factor * zeroCurl};
    }
    if (c < 0) {
        return {};
    }
    const Vec3 bcValue = edgeBasisValue(b, c, lambda, grad);
    const double bcCurl = 2.0 * dot(normal, cross(grad[static_cast<std::size_t>(b)], grad[static_cast<std::size_t>(c)]));
    if (ref.kind == LocalDofKind::FaceFirst0) {
        const double factor = lambda[static_cast<std::size_t>(c)];
        const Vec3 factorGrad = grad[static_cast<std::size_t>(c)];
        return {factor * zeroValue, dot(normal, cross(factorGrad, zeroValue)) + factor * zeroCurl};
    }
    // FaceFirst1: lambda_a * N_bc. Linearly independent from FaceFirst0
    // (lambda_c * N_ab); the previous form lambda_a*N_bc - lambda_b*N_ac
    // collapses algebraically to -lambda_c * N_ab = -FaceFirst0, which made
    // the port-face mass matrix rank-deficient and broke LAPACKE_dsygv on
    // basis-order 1 runs.
    const double factor = lambda[static_cast<std::size_t>(a)];
    const Vec3 factorGrad = grad[static_cast<std::size_t>(a)];
    return {factor * bcValue, dot(normal, cross(factorGrad, bcValue)) + factor * bcCurl};
}

}  // namespace

double poyntingPowerIntegral(const PortMode& mode, double frequencyHz) {
    const double omega = 2.0 * pi * frequencyHz;
    const double k0 = omega / c0;
    const double betaSq = k0 * k0 - mode.cutoffWavenumberSquared;
    if (betaSq <= 0.0 || mode.quadrature.empty()) {
        return 0.0;
    }
    const double beta = std::sqrt(betaSq);
    const double hScale = beta / (omega * mu0);
    double power = 0.0;
    for (const auto& qp : mode.quadrature) {
        const Vec3 E = qp.modeFieldValue;
        const Vec3 H = hScale * cross(qp.normal, E);
        power += 0.5 * qp.weight * dot(cross(E, H), qp.normal);
    }
    return power;
}

double powerNormalizationFactor(const PortMode& mode, double frequencyHz) {
    const double power = poyntingPowerIntegral(mode, frequencyHz);
    if (power <= 0.0) {
        return 0.0;
    }
    return 1.0 / std::sqrt(power);
}

PortModeSolver::PortModeSolver(const Mesh& mesh, const EdgeTopology& topology)
    : mesh_(mesh), topology_(topology) {
    const auto& edges = topology_.edges();
    for (std::size_t i = 0; i < edges.size(); ++i) {
        edgeDofsByNodePair_[edgeKey(edges[i].v0, edges[i].v1)] = {edges[i].zeroDof, edges[i].firstDof};
    }
}

const PortMode& PortModeSolver::solve(int faceId) const {
    const auto it = cache_.find(faceId);
    if (it != cache_.end()) {
        return it->second;
    }
    auto inserted = cache_.emplace(faceId, computeMode(faceId));
    return inserted.first->second;
}

void PortModeSolver::setPrecomputed(int faceId, PortMode mode) const {
    cache_[faceId] = std::move(mode);
}

void PortModeSolver::setMultiMode(int faceId, MultiPortMode multi) const {
    if (multi.modes.empty()) {
        throw std::runtime_error("PortModeSolver::setMultiMode: empty mode list for face " + std::to_string(faceId));
    }
    const int excIdx = std::clamp(multi.excitationModeIndex, 0, static_cast<int>(multi.modes.size()) - 1);
    multi.excitationModeIndex = excIdx;
    cache_[faceId] = multi.modes[static_cast<std::size_t>(excIdx)];
    multiCache_[faceId] = std::move(multi);
}

const MultiPortMode* PortModeSolver::multiMode(int faceId) const {
    const auto it = multiCache_.find(faceId);
    return (it == multiCache_.end()) ? nullptr : &it->second;
}

PortMode PortModeSolver::computeMode(int faceId) const {
    MultiPortMode multi = computeMultiMode(faceId, 1);
    if (multi.modes.empty()) {
        PortMode empty;
        empty.faceId = faceId;
        return empty;
    }
    return multi.modes.front();
}

MultiPortMode PortModeSolver::computeMultiMode(int faceId, int count) const {
    MultiPortMode multi;
    multi.faceId = faceId;
    multi.excitationModeIndex = 0;
    if (count < 1) {
        count = 1;
    }

    std::vector<SurfaceTriangle> triangles;
    for (const auto& tri : mesh_.surfaceTriangles) {
        const auto fit = mesh_.facets.find(tri.facetId);
        if (fit == mesh_.facets.end() || fit->second.faceId != faceId) {
            continue;
        }
        triangles.push_back(tri);
    }

    if (triangles.empty()) {
        return multi;
    }

    std::unordered_map<std::uint64_t, int> portEdgeUseCounts;
    for (const auto& tri : triangles) {
        const std::array<std::array<int, 2>, 3> localEdges{{{0, 1}, {1, 2}, {2, 0}}};
        for (const auto& localPair : localEdges) {
            const int aNode = tri.vertexIds[static_cast<std::size_t>(localPair[0])];
            const int bNode = tri.vertexIds[static_cast<std::size_t>(localPair[1])];
            ++portEdgeUseCounts[edgeKey(aNode, bNode)];
        }
    }

    std::map<int, int> portIndexByGlobal;
    std::vector<int> globalByPortIndex;
    std::vector<std::vector<PortDofRef>> triangleDofs;
    triangleDofs.reserve(triangles.size());
    auto portIndex = [&](int globalIndex) {
        const auto it = portIndexByGlobal.find(globalIndex);
        if (it != portIndexByGlobal.end()) {
            return it->second;
        }
        const int index = static_cast<int>(globalByPortIndex.size());
        portIndexByGlobal[globalIndex] = index;
        globalByPortIndex.push_back(globalIndex);
        return index;
    };

    for (const auto& tri : triangles) {
        std::vector<PortDofRef> refs;
        const std::array<std::array<int, 2>, 3> localEdges{{{0, 1}, {1, 2}, {2, 0}}};
        for (const auto& localPair : localEdges) {
            const int a = localPair[0];
            const int b = localPair[1];
            const int aNode = tri.vertexIds[static_cast<std::size_t>(a)];
            const int bNode = tri.vertexIds[static_cast<std::size_t>(b)];
            if (portEdgeUseCounts[edgeKey(aNode, bNode)] == 1) {
                continue;
            }
            const auto [zeroDof, firstDof] = edgeDofsForNodes(aNode, bNode);
            const int orientedA = aNode < bNode ? a : b;
            const int orientedB = aNode < bNode ? b : a;
            if (zeroDof >= 0) {
                refs.push_back({portIndex(zeroDof), zeroDof, LocalDofKind::EdgeZero, {{orientedA, orientedB, -1}}, 1.0});
            }
            if (topology_.basisOrder() >= 1 && firstDof >= 0) {
                refs.push_back({portIndex(firstDof), firstDof, LocalDofKind::EdgeFirst, {{orientedA, orientedB, -1}}, 1.0});
            }
        }
        if (topology_.basisOrder() >= 1) {
            std::array<std::pair<int, int>, 3> sortedNodes{{{tri.vertexIds[0], 0}, {tri.vertexIds[1], 1}, {tri.vertexIds[2], 2}}};
            std::sort(sortedNodes.begin(), sortedNodes.end(), [](const auto& lhs, const auto& rhs) {
                return lhs.first < rhs.first;
            });
            const std::array<int, 3> sortedLocal{{sortedNodes[0].second, sortedNodes[1].second, sortedNodes[2].second}};
            const auto faceDofs = topology_.faceFirstDofsForNodes(tri.vertexIds[0], tri.vertexIds[1], tri.vertexIds[2]);
            for (int i = 0; i < 2; ++i) {
                const int dof = faceDofs[static_cast<std::size_t>(i)];
                if (dof >= 0) {
                    refs.push_back({portIndex(dof), dof, i == 0 ? LocalDofKind::FaceFirst0 : LocalDofKind::FaceFirst1, sortedLocal, 1.0});
                }
            }
        }
        triangleDofs.push_back(std::move(refs));
    }

    const std::size_t n = globalByPortIndex.size();
    if (n == 0) {
        return multi;
    }

    Matrix stiffness = makeMatrix(n);
    Matrix mass = makeMatrix(n);
    const auto quadrature = triangleQuadraturePoints();
    // Cache per-triangle (p, area, normal, grad) for reuse in the eigensolve
    // and for quadrature-point reconstruction afterwards.
    struct TriCache {
        std::array<Vec3, 3> p{};
        double area = 0.0;
        Vec3 normal{};
        std::array<Vec3, 3> grad{};
    };
    std::vector<TriCache> triCache(triangles.size());
    for (std::size_t triIndex = 0; triIndex < triangles.size(); ++triIndex) {
        const auto& tri = triangles[triIndex];
        std::array<Vec3, 3> p{};
        for (int i = 0; i < 3; ++i) {
            p[static_cast<std::size_t>(i)] = mesh_.pointsById.at(tri.vertexIds[static_cast<std::size_t>(i)]);
        }
        double area = 0.0;
        Vec3 normal;
        const auto grad = surfaceGradients(p, area, normal);
        triCache[triIndex] = TriCache{p, area, normal, grad};
        if (area <= 1.0e-30) {
            continue;
        }
        const auto& refs = triangleDofs[triIndex];
        for (std::size_t a = 0; a < refs.size(); ++a) {
            for (std::size_t b = 0; b < refs.size(); ++b) {
                double localStiffness = 0.0;
                double localMass = 0.0;
                for (const auto& qp : quadrature) {
                    const BasisValue rowBasis = evaluateBasis(refs[a], qp.lambda, grad, normal);
                    const BasisValue colBasis = evaluateBasis(refs[b], qp.lambda, grad, normal);
                    localStiffness += area * qp.weight * rowBasis.curl * colBasis.curl;
                    localMass += area * qp.weight * dot(rowBasis.value, colBasis.value);
                }
                stiffness[static_cast<std::size_t>(refs[a].portIndex)][static_cast<std::size_t>(refs[b].portIndex)] += localStiffness;
                mass[static_cast<std::size_t>(refs[a].portIndex)][static_cast<std::size_t>(refs[b].portIndex)] += localMass;
            }
        }
    }

    const auto [eigenvalues, eigenvectors] = generalizedEigen(stiffness, mass);

    // Sort eigenvalues ascending. The H(curl) generalized eigenproblem
    // K v = lambda M v on a closed port has a non-trivial null space coming
    // from the gradient subspace (curl ker), giving spurious near-zero
    // eigenvalues that are not physical TE/TM modes. We filter them out and
    // keep the lowest `count` strictly-positive eigenvalues.
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return eigenvalues[a] < eigenvalues[b]; });
    const double maxEigenvalue = std::abs(eigenvalues[order.back()]);
    const double threshold = std::max(1.0e-9, maxEigenvalue * 1.0e-10);

    std::vector<std::size_t> selectedIndices;
    selectedIndices.reserve(static_cast<std::size_t>(count));
    for (std::size_t idx : order) {
        if (eigenvalues[idx] > threshold) {
            selectedIndices.push_back(idx);
            if (static_cast<int>(selectedIndices.size()) >= count) {
                break;
            }
        }
    }
    if (selectedIndices.empty()) {
        // Pathological: keep at least the dominant entry so downstream code
        // gets a usable PortMode (matching legacy behavior).
        selectedIndices.push_back(order.front());
    }

    multi.modes.reserve(selectedIndices.size());
    for (std::size_t selected : selectedIndices) {
        std::vector<double> dofs(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            dofs[i] = eigenvectors[i][selected];
        }
        auto massDofs = multiply(mass, dofs);
        double normM = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            normM += dofs[i] * massDofs[i];
        }
        normM = std::sqrt(std::max(normM, 1.0e-30));
        for (double& value : dofs) {
            value /= normM;
        }
        massDofs = multiply(mass, dofs);

        PortMode pm;
        pm.faceId = faceId;
        pm.cutoffWavenumberSquared = std::max(0.0, eigenvalues[selected]);
        for (std::size_t i = 0; i < dofs.size(); ++i) {
            if (std::abs(dofs[i]) > 1.0e-18) {
                pm.edgeDofs.emplace_back(globalByPortIndex[i], dofs[i]);
            }
            if (std::abs(massDofs[i]) > 1.0e-18) {
                pm.couplingWeights.emplace_back(globalByPortIndex[i], massDofs[i]);
            }
        }
        for (std::size_t triIndex = 0; triIndex < triangles.size(); ++triIndex) {
            const auto& cache = triCache[triIndex];
            if (cache.area <= 1.0e-30) {
                continue;
            }
            const auto& refs = triangleDofs[triIndex];
            for (const auto& qp : quadrature) {
                Vec3 fieldValue{};
                for (const auto& ref : refs) {
                    const BasisValue basis = evaluateBasis(ref, qp.lambda, cache.grad, cache.normal);
                    fieldValue = fieldValue + dofs[static_cast<std::size_t>(ref.portIndex)] * basis.value;
                }
                pm.quadrature.push_back(PortQuadraturePoint{fieldValue, cache.normal, cache.area * qp.weight});
            }
        }
        multi.modes.push_back(std::move(pm));
    }

    return multi;
}

std::pair<int, int> PortModeSolver::edgeDofsForNodes(int aNode, int bNode) const {
    const auto it = edgeDofsByNodePair_.find(edgeKey(aNode, bNode));
    if (it == edgeDofsByNodePair_.end()) {
        return {-1, -1};
    }
    return it->second;
}

std::uint64_t PortModeSolver::edgeKey(int aNode, int bNode) {
    const int lo = std::min(aNode, bNode);
    const int hi = std::max(aNode, bNode);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo)) << 32U) | static_cast<std::uint32_t>(hi);
}

}  // namespace fem
