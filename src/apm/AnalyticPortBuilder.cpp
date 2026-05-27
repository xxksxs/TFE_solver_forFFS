#include "bpfem/apm/AnalyticPortBuilder.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/core/Math.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fem::apm {

namespace {

// Wandzura 7-point rule on the reference triangle (degree 5). Identical to the
// rule used by PortModeSolver so projections of the analytic mode reach the
// same accuracy as the numerical eigenmode.
struct TriangleQuadraturePoint {
    std::array<double, 3> lambda;
    double weight = 0.0;  // sum of weights = 1
};

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

std::array<Vec3, 3> surfaceGradients(const std::array<Vec3, 3>& p, double& area, Vec3& normal) {
    const Vec3 e01 = p[1] - p[0];
    const Vec3 e02 = p[2] - p[0];
    Vec3 raw = cross(e01, e02);
    const double len = norm(raw);
    area = 0.5 * len;
    if (len <= 1.0e-30 || norm(e01) <= 1.0e-30) {
        return {};
    }
    normal = raw / len;
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

struct PortDofRef {
    int globalIndex = -1;
    LocalDofKind kind = LocalDofKind::EdgeZero;
    std::array<int, 3> localNodes{{-1, -1, -1}};
};

Vec3 evaluatePortBasis(const PortDofRef& ref,
                      const std::array<double, 3>& lambda,
                      const std::array<Vec3, 3>& grad) {
    const int a = ref.localNodes[0];
    const int b = ref.localNodes[1];
    const int c = ref.localNodes[2];
    if (a < 0 || b < 0) {
        return {};
    }
    const Vec3 zeroValue = edgeBasisValue(a, b, lambda, grad);
    if (ref.kind == LocalDofKind::EdgeZero) {
        return zeroValue;
    }
    if (ref.kind == LocalDofKind::EdgeFirst) {
        const double factor = lambda[static_cast<std::size_t>(a)] - lambda[static_cast<std::size_t>(b)];
        return factor * zeroValue;
    }
    if (c < 0) {
        return {};
    }
    if (ref.kind == LocalDofKind::FaceFirst0) {
        return lambda[static_cast<std::size_t>(c)] * zeroValue;
    }
    const Vec3 bcValue = edgeBasisValue(b, c, lambda, grad);
    return lambda[static_cast<std::size_t>(a)] * bcValue;
}

std::uint64_t edgeKey(int aNode, int bNode) {
    const int lo = std::min(aNode, bNode);
    const int hi = std::max(aNode, bNode);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo)) << 32U) | static_cast<std::uint32_t>(hi);
}

std::pair<double, double> toUV(const RectPortGeometry& g, const Vec3& p) {
    const Vec3 d = p - g.origin;
    return {dot(d, g.axisU), dot(d, g.axisV)};
}

// TE10 analytic mode in port-local UV coordinates with a >= b:
//   e_t(u, v) = sqrt(2 / (a*b)) * sin(pi * u / a) * v_hat
Vec3 analyticTE10(const RectPortGeometry& g, double u, double v) {
    (void)v;
    const double amp = std::sqrt(2.0 / (g.a * g.b)) * std::sin(pi * u / g.a);
    return amp * g.axisV;
}

}  // namespace

AnalyticPortBuilder::AnalyticPortBuilder(const Mesh& mesh, const EdgeTopology& topology)
    : mesh_(mesh), topology_(topology) {}

RectPortGeometry AnalyticPortBuilder::detectRectangular(int faceId) const {
    RectPortGeometry geo;
    std::vector<SurfaceTriangle> triangles;
    for (const auto& tri : mesh_.surfaceTriangles) {
        const auto fit = mesh_.facets.find(tri.facetId);
        if (fit == mesh_.facets.end() || fit->second.faceId != faceId) {
            continue;
        }
        triangles.push_back(tri);
    }
    if (triangles.empty()) {
        return geo;
    }

    Vec3 normalAccum{0.0, 0.0, 0.0};
    Vec3 centroidAccum{0.0, 0.0, 0.0};
    double totalArea = 0.0;
    std::map<int, Vec3> uniqueVertices;
    for (const auto& tri : triangles) {
        std::array<Vec3, 3> p{};
        for (int i = 0; i < 3; ++i) {
            const int vid = tri.vertexIds[static_cast<std::size_t>(i)];
            const auto pit = mesh_.pointsById.find(vid);
            if (pit == mesh_.pointsById.end()) {
                return geo;
            }
            p[static_cast<std::size_t>(i)] = pit->second;
            uniqueVertices.emplace(vid, pit->second);
        }
        const Vec3 e01 = p[1] - p[0];
        const Vec3 e02 = p[2] - p[0];
        const Vec3 n = cross(e01, e02);
        const double area = 0.5 * norm(n);
        normalAccum = normalAccum + n;
        centroidAccum = centroidAccum + area * (p[0] + p[1] + p[2]) / 3.0;
        totalArea += area;
    }
    if (totalArea <= 1.0e-30) {
        return geo;
    }
    const double normalLength = norm(normalAccum);
    if (normalLength <= 1.0e-30) {
        return geo;
    }
    geo.normal = normalAccum / normalLength;
    Vec3 centroid = centroidAccum / totalArea;

    std::vector<Vec3> projected;
    projected.reserve(uniqueVertices.size());
    for (const auto& [vid, p] : uniqueVertices) {
        const Vec3 r = p - centroid;
        const Vec3 r_par = r - dot(r, geo.normal) * geo.normal;
        projected.push_back(r_par);
    }
    if (projected.size() < 3) {
        return geo;
    }

    auto orthBasis = [&]() {
        Vec3 ref = (std::abs(geo.normal.x) < 0.9) ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
        Vec3 u0 = ref - dot(ref, geo.normal) * geo.normal;
        u0 = u0 / std::max(norm(u0), 1.0e-30);
        const Vec3 v0 = cross(geo.normal, u0);
        return std::pair<Vec3, Vec3>{u0, v0};
    };
    auto [u0, v0] = orthBasis();
    double cuu = 0.0, cvv = 0.0, cuv = 0.0;
    for (const auto& r : projected) {
        const double pu = dot(r, u0);
        const double pv = dot(r, v0);
        cuu += pu * pu;
        cvv += pv * pv;
        cuv += pu * pv;
    }
    const double phi = 0.5 * std::atan2(2.0 * cuv, cuu - cvv);
    const Vec3 axis1 = std::cos(phi) * u0 + std::sin(phi) * v0;
    const Vec3 axis2 = -std::sin(phi) * u0 + std::cos(phi) * v0;

    double minA = std::numeric_limits<double>::max();
    double maxA = std::numeric_limits<double>::lowest();
    double minB = std::numeric_limits<double>::max();
    double maxB = std::numeric_limits<double>::lowest();
    for (const auto& r : projected) {
        const double pa = dot(r, axis1);
        const double pb = dot(r, axis2);
        minA = std::min(minA, pa);
        maxA = std::max(maxA, pa);
        minB = std::min(minB, pb);
        maxB = std::max(maxB, pb);
    }
    const double extentA = maxA - minA;
    const double extentB = maxB - minB;
    if (extentA <= 0.0 || extentB <= 0.0) {
        return geo;
    }
    if (extentA >= extentB) {
        geo.axisU = axis1;
        geo.axisV = axis2;
        geo.a = extentA;
        geo.b = extentB;
        geo.origin = centroid + minA * axis1 + minB * axis2;
    } else {
        geo.axisU = axis2;
        geo.axisV = axis1;
        geo.a = extentB;
        geo.b = extentA;
        geo.origin = centroid + minB * axis2 + minA * axis1;
    }
    geo.valid = true;
    return geo;
}

PortMode AnalyticPortBuilder::build(int faceId) const {
    PortMode mode;
    mode.faceId = faceId;

    const RectPortGeometry geo = detectRectangular(faceId);
    if (!geo.valid) {
        throw std::runtime_error("AnalyticPortBuilder: failed to detect rectangular geometry on face " + std::to_string(faceId));
    }
    mode.cutoffWavenumberSquared = (pi / geo.a) * (pi / geo.a);

    std::vector<SurfaceTriangle> triangles;
    for (const auto& tri : mesh_.surfaceTriangles) {
        const auto fit = mesh_.facets.find(tri.facetId);
        if (fit == mesh_.facets.end() || fit->second.faceId != faceId) {
            continue;
        }
        triangles.push_back(tri);
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
    auto portIndexOf = [&](int globalIndex) {
        const auto it = portIndexByGlobal.find(globalIndex);
        if (it != portIndexByGlobal.end()) {
            return it->second;
        }
        const int idx = static_cast<int>(globalByPortIndex.size());
        portIndexByGlobal[globalIndex] = idx;
        globalByPortIndex.push_back(globalIndex);
        return idx;
    };
    struct LocalRef {
        PortDofRef ref;
        int portIndex = -1;
    };
    std::vector<std::vector<LocalRef>> triangleLocalRefs;
    triangleLocalRefs.reserve(triangles.size());

    for (const auto& tri : triangles) {
        std::vector<LocalRef> refs;
        const std::array<std::array<int, 2>, 3> localEdges{{{0, 1}, {1, 2}, {2, 0}}};
        for (const auto& localPair : localEdges) {
            const int aLocal = localPair[0];
            const int bLocal = localPair[1];
            const int aNode = tri.vertexIds[static_cast<std::size_t>(aLocal)];
            const int bNode = tri.vertexIds[static_cast<std::size_t>(bLocal)];
            if (portEdgeUseCounts[edgeKey(aNode, bNode)] == 1) {
                continue;
            }
            const auto [zeroDof, firstDof] = topology_.edgeDofsForNodes(aNode, bNode);
            const int orientedA = aNode < bNode ? aLocal : bLocal;
            const int orientedB = aNode < bNode ? bLocal : aLocal;
            if (zeroDof >= 0) {
                refs.push_back({{zeroDof, LocalDofKind::EdgeZero, {{orientedA, orientedB, -1}}}, portIndexOf(zeroDof)});
            }
            if (topology_.basisOrder() >= 1 && firstDof >= 0) {
                refs.push_back({{firstDof, LocalDofKind::EdgeFirst, {{orientedA, orientedB, -1}}}, portIndexOf(firstDof)});
            }
        }
        if (topology_.basisOrder() >= 1) {
            std::array<std::pair<int, int>, 3> sortedNodes{{{tri.vertexIds[0], 0}, {tri.vertexIds[1], 1}, {tri.vertexIds[2], 2}}};
            std::sort(sortedNodes.begin(), sortedNodes.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
            const std::array<int, 3> sortedLocal{{sortedNodes[0].second, sortedNodes[1].second, sortedNodes[2].second}};
            const auto faceDofs = topology_.faceFirstDofsForNodes(tri.vertexIds[0], tri.vertexIds[1], tri.vertexIds[2]);
            for (int i = 0; i < 2; ++i) {
                const int dof = faceDofs[static_cast<std::size_t>(i)];
                if (dof >= 0) {
                    refs.push_back({{dof, i == 0 ? LocalDofKind::FaceFirst0 : LocalDofKind::FaceFirst1, sortedLocal}, portIndexOf(dof)});
                }
            }
        }
        triangleLocalRefs.push_back(std::move(refs));
    }

    const int n = static_cast<int>(globalByPortIndex.size());
    if (n == 0) {
        return mode;
    }

    std::vector<std::vector<double>> M(static_cast<std::size_t>(n), std::vector<double>(static_cast<std::size_t>(n), 0.0));
    std::vector<double> b(static_cast<std::size_t>(n), 0.0);

    const auto quadrature = triangleQuadraturePoints();
    for (std::size_t tIndex = 0; tIndex < triangles.size(); ++tIndex) {
        const auto& tri = triangles[tIndex];
        std::array<Vec3, 3> p{};
        for (int i = 0; i < 3; ++i) {
            const auto pit = mesh_.pointsById.find(tri.vertexIds[static_cast<std::size_t>(i)]);
            if (pit == mesh_.pointsById.end()) {
                continue;
            }
            p[static_cast<std::size_t>(i)] = pit->second;
        }
        double area = 0.0;
        Vec3 normal{};
        const auto grad = surfaceGradients(p, area, normal);
        if (area <= 1.0e-30) {
            continue;
        }
        const auto& refs = triangleLocalRefs[tIndex];

        for (const auto& qp : quadrature) {
            const double w = area * qp.weight;
            const Vec3 xyz = qp.lambda[0] * p[0] + qp.lambda[1] * p[1] + qp.lambda[2] * p[2];
            const auto [u, v] = toUV(geo, xyz);
            const Vec3 et = analyticTE10(geo, u, v);
            std::vector<Vec3> N(refs.size());
            for (std::size_t i = 0; i < refs.size(); ++i) {
                N[i] = evaluatePortBasis(refs[i].ref, qp.lambda, grad);
            }
            for (std::size_t i = 0; i < refs.size(); ++i) {
                b[static_cast<std::size_t>(refs[i].portIndex)] += w * dot(et, N[i]);
                for (std::size_t j = 0; j < refs.size(); ++j) {
                    M[static_cast<std::size_t>(refs[i].portIndex)][static_cast<std::size_t>(refs[j].portIndex)] += w * dot(N[i], N[j]);
                }
            }
        }
    }

    {
        std::vector<std::vector<double>> L(static_cast<std::size_t>(n), std::vector<double>(static_cast<std::size_t>(n), 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = M[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                for (int k = 0; k < j; ++k) {
                    sum -= L[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] * L[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)];
                }
                if (i == j) {
                    if (sum < 1.0e-30) {
                        throw std::runtime_error("AnalyticPortBuilder: port-face mass matrix not SPD on face " + std::to_string(faceId));
                    }
                    L[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = std::sqrt(sum);
                } else {
                    L[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = sum / L[static_cast<std::size_t>(j)][static_cast<std::size_t>(j)];
                }
            }
        }
        std::vector<double> y(static_cast<std::size_t>(n), 0.0);
        for (int i = 0; i < n; ++i) {
            double sum = b[static_cast<std::size_t>(i)];
            for (int k = 0; k < i; ++k) {
                sum -= L[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] * y[static_cast<std::size_t>(k)];
            }
            y[static_cast<std::size_t>(i)] = sum / L[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
        }
        std::vector<double> dofs(static_cast<std::size_t>(n), 0.0);
        for (int i = n - 1; i >= 0; --i) {
            double sum = y[static_cast<std::size_t>(i)];
            for (int k = i + 1; k < n; ++k) {
                sum -= L[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] * dofs[static_cast<std::size_t>(k)];
            }
            dofs[static_cast<std::size_t>(i)] = sum / L[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
        }

        double l2sq = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                l2sq += dofs[static_cast<std::size_t>(i)] * M[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] * dofs[static_cast<std::size_t>(j)];
            }
        }
        const double alpha = std::sqrt(std::max(l2sq, 1.0e-300));
        for (auto& d : dofs) {
            d /= alpha;
        }
        for (auto& bi : b) {
            bi /= alpha;
        }

        for (int i = 0; i < n; ++i) {
            const double w = b[static_cast<std::size_t>(i)];
            if (std::abs(w) > 1.0e-18) {
                mode.couplingWeights.emplace_back(globalByPortIndex[static_cast<std::size_t>(i)], w);
            }
            const double d = dofs[static_cast<std::size_t>(i)];
            if (std::abs(d) > 1.0e-18) {
                mode.edgeDofs.emplace_back(globalByPortIndex[static_cast<std::size_t>(i)], d);
            }
        }

        for (std::size_t tIndex = 0; tIndex < triangles.size(); ++tIndex) {
            const auto& tri = triangles[tIndex];
            std::array<Vec3, 3> p{};
            for (int i = 0; i < 3; ++i) {
                const auto pit = mesh_.pointsById.find(tri.vertexIds[static_cast<std::size_t>(i)]);
                if (pit == mesh_.pointsById.end()) {
                    continue;
                }
                p[static_cast<std::size_t>(i)] = pit->second;
            }
            double area = 0.0;
            Vec3 normal{};
            const auto grad = surfaceGradients(p, area, normal);
            if (area <= 1.0e-30) {
                continue;
            }
            const auto& refs = triangleLocalRefs[tIndex];
            for (const auto& qp : quadrature) {
                Vec3 fieldValue{0.0, 0.0, 0.0};
                for (const auto& lr : refs) {
                    const Vec3 N = evaluatePortBasis(lr.ref, qp.lambda, grad);
                    fieldValue = fieldValue + dofs[static_cast<std::size_t>(lr.portIndex)] * N;
                }
                mode.quadrature.push_back(PortQuadraturePoint{fieldValue, normal, area * qp.weight});
            }
        }
    }

    return mode;
}

}  // namespace fem::apm
