#include "bpfem/post/OutputWriter.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/core/Math.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace fem {

namespace {

struct BasisValue {
    Vec3 value;
    Vec3 curl;
};

// H(curl) hierarchical basis on a tetrahedron, evaluated at barycentric
// coordinates lambda. Returns both the vector value (for E reconstruction) and
// the curl (for H = j curl(E) / (omega mu_r mu_0)). The forms here must match
// FEMAssembler::evaluateBasis exactly so that the visualization is consistent
// with the assembled system.
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
        // FaceFirst0 = lambda_c * N_ab
        const double factor = lambda[static_cast<std::size_t>(c)];
        const Vec3 factorGrad = grad[static_cast<std::size_t>(c)];
        return {factor * zeroValue, cross(factorGrad, zeroValue) + factor * zeroCurl};
    }
    // FaceFirst1 = lambda_a * N_bc (linearly independent from FaceFirst0; see
    // PortModeSolver / FEMAssembler for the matching definition).
    const double factor = lambda[static_cast<std::size_t>(a)];
    const Vec3 factorGrad = grad[static_cast<std::size_t>(a)];
    return {factor * bcValue, cross(factorGrad, bcValue) + factor * bcCurl};
}

// VTK Lagrange tetrahedron node ordering (cell type 71, vtkLagrangeTetra).
// The same corner+edge ordering also matches VTK_QUADRATIC_TETRA (cell type
// 24) at order 2, so the table built below is reusable for both writers.
// Edge order (matches vtkTetra::edges_ and VTK_QUADRATIC_TETRA midpoints):
//   E0=(0,1) E1=(1,2) E2=(2,0) E3=(0,3) E4=(1,3) E5=(2,3)
// Face order (matches vtkTetra::faces_):
//   F0=(0,1,3) F1=(1,2,3) F2=(2,0,3) F3=(0,2,1)
struct EdgeKey { int a; int b; };
struct FaceKey { int a; int b; int c; };

constexpr std::array<EdgeKey, 6> kVtkEdges = {{
    {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}
}};
constexpr std::array<FaceKey, 4> kVtkFaces = {{
    {0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 2, 1}
}};

int classifyOnFace(const std::array<double, 4>& lambda) {
    constexpr double tol = 1.0e-12;
    int zeroIdx = -1;
    for (int k = 0; k < 4; ++k) {
        if (lambda[static_cast<std::size_t>(k)] < tol) {
            if (zeroIdx >= 0) {
                return -1;  // multiple zeros => corner or edge interior
            }
            zeroIdx = k;
        }
    }
    return zeroIdx;
}

Material materialForBody(const Mesh& mesh, const ProjectDefinition& project, int bodyId) {
    auto bit = mesh.bodies.find(bodyId);
    if (bit != mesh.bodies.end()) {
        auto mit = project.materials.find(bit->second.name);
        if (mit != project.materials.end()) {
            return mit->second;
        }
    }
    auto mit = project.materials.find(project.backgroundMaterial);
    if (mit != project.materials.end()) {
        return mit->second;
    }
    return Material{"vacuum", 1.0, 1.0, 0.0};
}

}  // namespace

OutputWriter::LagrangeTetTable OutputWriter::buildLagrangeTetTable(int order) {
    if (order < kMinFieldOutputOrder || order > kMaxFieldOutputOrder) {
        throw std::runtime_error("buildLagrangeTetTable: order out of range");
    }
    LagrangeTetTable table;
    table.order = order;

    auto append = [&](const std::array<double, 4>& lambda,
                      int faceIndex,
                      int cornerIndex,
                      int localEdge,
                      int localFaceCenter) {
        table.lambda.push_back(lambda);
        table.faceOnly.push_back(faceIndex);
        table.cornerIndex.push_back(cornerIndex);
        table.localEdge.push_back(localEdge);
        table.localFaceCenter.push_back(localFaceCenter);
    };

    // Corners (in VTK order: 0,1,2,3).
    for (int corner = 0; corner < 4; ++corner) {
        std::array<double, 4> lambda{{0.0, 0.0, 0.0, 0.0}};
        lambda[static_cast<std::size_t>(corner)] = 1.0;
        append(lambda, -1, corner, -1, -1);
    }

    if (order == 1) {
        return table;
    }

    // Edge interior nodes: for each VTK edge (a,b) emit (p-1) equispaced points
    // strictly between corner a and corner b. At order 2 there is one edge
    // midpoint per edge, so localEdge maps each emitted node to its edge index.
    // At order 3 each edge gets 2 interior points; both still belong to the
    // same edge so they share the edge id and the per-cell midpoint sharing
    // covers the central one. (At order 3 the two interior points are not at
    // the same xyz on the shared edge across cells, so VTK convention treats
    // each as its own node; for sharing we still group by edge index because
    // the order-3 node t-coordinate is the same on either side of the edge.)
    for (std::size_t e = 0; e < kVtkEdges.size(); ++e) {
        const auto& edge = kVtkEdges[e];
        for (int k = 1; k <= order - 1; ++k) {
            const double t = static_cast<double>(k) / static_cast<double>(order);
            std::array<double, 4> lambda{{0.0, 0.0, 0.0, 0.0}};
            lambda[static_cast<std::size_t>(edge.a)] = 1.0 - t;
            lambda[static_cast<std::size_t>(edge.b)] = t;
            append(lambda, classifyOnFace(lambda), -1, static_cast<int>(e), -1);
        }
    }

    if (order < 3) {
        return table;
    }

    // Face interior nodes: barycentric (i, j, k) with i+j+k = p, all > 0.
    // At order 3 each face has exactly one interior node (i=j=k=1, scaled by
    // 1/3), shared between the (up to) two cells incident on that face.
    for (std::size_t f = 0; f < kVtkFaces.size(); ++f) {
        const auto& face = kVtkFaces[f];
        for (int j = 1; j <= order - 2; ++j) {
            for (int i = 1; i <= order - 1 - j; ++i) {
                const int k = order - i - j;
                if (k < 1) {
                    continue;
                }
                std::array<double, 4> lambda{{0.0, 0.0, 0.0, 0.0}};
                lambda[static_cast<std::size_t>(face.a)] = static_cast<double>(i) / order;
                lambda[static_cast<std::size_t>(face.b)] = static_cast<double>(j) / order;
                lambda[static_cast<std::size_t>(face.c)] = static_cast<double>(k) / order;
                append(lambda, classifyOnFace(lambda), -1, -1, static_cast<int>(f));
            }
        }
    }

    // Volume interior nodes: barycentric (i, j, k, l) with i+j+k+l = p, all > 0.
    for (int l = 1; l <= order - 3; ++l) {
        for (int k = 1; k <= order - 2 - l + 1; ++k) {
            for (int j = 1; j <= order - 1 - l - k + 1; ++j) {
                const int i = order - j - k - l;
                if (i < 1) {
                    continue;
                }
                std::array<double, 4> lambda{
                    static_cast<double>(i) / order,
                    static_cast<double>(j) / order,
                    static_cast<double>(k) / order,
                    static_cast<double>(l) / order};
                append(lambda, -1, -1, -1, -1);
            }
        }
    }

    return table;
}

void OutputWriter::writeSParameters(const std::filesystem::path& path, const std::vector<SParameterPoint>& data) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot write S-parameter file: " + path.string());
    }
    auto db = [](const std::complex<double>& value) {
        return 20.0 * std::log10(std::max(std::abs(value), 1.0e-300));
    };
    out << "freq_Hz,S11_real,S11_imag,S11_dB,S21_real,S21_imag,S21_dB\n";
    out << std::setprecision(15);
    for (const auto& row : data) {
        out << row.frequencyHz << ','
            << row.s11.real() << ',' << row.s11.imag() << ',' << db(row.s11) << ','
            << row.s21.real() << ',' << row.s21.imag() << ',' << db(row.s21) << '\n';
    }
}

void OutputWriter::writeVTU(const std::filesystem::path& path,
                            const Mesh& mesh,
                            const ProjectDefinition& project,
                            const EdgeTopology& topology,
                            const std::vector<std::complex<double>>& edgeDofs,
                            double frequencyHz,
                            int fieldOutputOrder) {
    if (fieldOutputOrder < kMinFieldOutputOrder || fieldOutputOrder > kMaxFieldOutputOrder) {
        throw std::runtime_error("writeVTU: fieldOutputOrder must be in [1, 3]");
    }
#ifndef BPFEM_OUTPUT_PARAVIEW
    if (fieldOutputOrder > 2) {
        throw std::runtime_error(
            "writeVTU: --field-output-order 3 requires the ParaView VTK Lagrange writer. "
            "Rebuild with -DBPFEM_OUTPUT_PARAVIEW=ON or pass --field-output-order 2.");
    }
#endif
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot write VTU file: " + path.string());
    }
    const auto table = buildLagrangeTetTable(fieldOutputOrder);
    const std::size_t nodesPerCell = table.lambda.size();
    const std::size_t totalCells = mesh.tetrahedra.size();

    // ----- Shared-node mapping -----
    //
    // Nedelec H(curl) elements are tangentially continuous across cell
    // boundaries by construction, so values reconstructed at *shared* Lagrange
    // sample positions (corners on cell-cell edges, edge midpoints on shared
    // edges, face-center points on shared faces) are physically equal on
    // either side up to the small normal-jump. We collapse them into a single
    // VTU node per geometric entity, which shrinks an order-2 file from
    // 10 nodes/cell down to ~ 1.7 nodes/cell on tet meshes.
    //
    // Mapping rules per Lagrange node n in a cell:
    //   - corner       (table.cornerIndex[n] >= 0)         -> shared by
    //                  vertex id;  global id = vertex slot.
    //   - edge interior (table.localEdge[n] >= 0)           -> shared by
    //                  geometric edge id (from EdgeTopology). At order 2
    //                  there is one midpoint per edge so the slot is unique;
    //                  at order 3 each edge has two interior points and we
    //                  encode them as (edge_id * 2 + k) where k = 0 or 1
    //                  along the lower-vertex-to-higher-vertex direction.
    //   - face interior (table.localFaceCenter[n] >= 0)     -> shared by
    //                  triangular face (sorted-vertex triple). Only at order 3.
    //   - volume interior                                  -> cell-local.
    auto edgeKey = [](int a, int b) -> std::uint64_t {
        const int lo = std::min(a, b);
        const int hi = std::max(a, b);
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo)) << 32U)
             | static_cast<std::uint32_t>(hi);
    };
    auto faceKey = [](int a, int b, int c) -> std::array<int, 3> {
        std::array<int, 3> v{a, b, c};
        std::sort(v.begin(), v.end());
        return v;
    };

    // Geometric-edge index map.
    std::unordered_map<std::uint64_t, int> edgeIndexByKey;
    edgeIndexByKey.reserve(topology.edges().size() * 2);
    {
        const auto& edges = topology.edges();
        for (std::size_t i = 0; i < edges.size(); ++i) {
            edgeIndexByKey[edgeKey(edges[i].v0, edges[i].v1)] = static_cast<int>(i);
        }
    }
    const int geometricEdgeCount = static_cast<int>(topology.edges().size());

    // Vertex-id -> contiguous slot. We use the pointId order as written by
    // the mesh; Mesh::pointIds is already sorted.
    std::unordered_map<int, int> vertexSlotById;
    vertexSlotById.reserve(mesh.pointIds.size() * 2);
    for (std::size_t i = 0; i < mesh.pointIds.size(); ++i) {
        vertexSlotById[mesh.pointIds[i]] = static_cast<int>(i);
    }

    // Allocate base ranges:
    //   [0, vertexCount)                                                 -> corners
    //   [vertexCount, vertexCount + edgeCount * (order-1))                -> edge points
    //   [..., + faceCount * faceInteriorPerFace)                          -> face centers (order=3)
    //   [..., totalGlobalNodes)                                          -> per-cell volume interiors
    const int vertexCount = static_cast<int>(mesh.pointIds.size());
    const int edgePointsPerEdge = std::max(0, fieldOutputOrder - 1);  // 0/1/2 for order 1/2/3
    const int edgeBase = vertexCount;
    const int edgeBlockSize = geometricEdgeCount * edgePointsPerEdge;
    const int faceBase = edgeBase + edgeBlockSize;

    // Face mapping: build by walking cells and sorting each face triple.
    // At order 3 each face has one interior point.
    const int faceInteriorPerFace = (fieldOutputOrder >= 3) ? 1 : 0;
    std::map<std::array<int, 3>, int> faceSlotByKey;  // unique face -> slot index
    if (faceInteriorPerFace > 0) {
        for (const auto& tet : mesh.tetrahedra) {
            for (const auto& f : kVtkFaces) {
                const std::array<int, 3> tri = faceKey(
                    tet.vertexIds[static_cast<std::size_t>(f.a)],
                    tet.vertexIds[static_cast<std::size_t>(f.b)],
                    tet.vertexIds[static_cast<std::size_t>(f.c)]);
                faceSlotByKey.emplace(tri, static_cast<int>(faceSlotByKey.size()));
            }
        }
    }
    const int faceBlockSize = static_cast<int>(faceSlotByKey.size()) * faceInteriorPerFace;
    int volumeBase = faceBase + faceBlockSize;

    // Count per-cell volume-interior nodes (those with all classification
    // fields == -1) and assign global ids.
    int volumeInteriorPerCell = 0;
    for (std::size_t n = 0; n < nodesPerCell; ++n) {
        if (table.cornerIndex[n] < 0
            && table.localEdge[n] < 0
            && table.localFaceCenter[n] < 0) {
            ++volumeInteriorPerCell;
        }
    }
    const int totalGlobalNodes = volumeBase + static_cast<int>(totalCells) * volumeInteriorPerCell;

    // Build cellLocalToGlobal map.
    std::vector<std::vector<int>> cellLocalToGlobal(totalCells, std::vector<int>(nodesPerCell, -1));
    {
        int volumeCursor = volumeBase;
        for (std::size_t cellIndex = 0; cellIndex < totalCells; ++cellIndex) {
            const auto& tet = mesh.tetrahedra[cellIndex];
            for (std::size_t n = 0; n < nodesPerCell; ++n) {
                int globalId = -1;
                if (table.cornerIndex[n] >= 0) {
                    const int vId = tet.vertexIds[static_cast<std::size_t>(table.cornerIndex[n])];
                    const auto vit = vertexSlotById.find(vId);
                    if (vit != vertexSlotById.end()) {
                        globalId = vit->second;
                    }
                } else if (table.localEdge[n] >= 0 && edgePointsPerEdge > 0) {
                    const auto& e = kVtkEdges[static_cast<std::size_t>(table.localEdge[n])];
                    const int va = tet.vertexIds[static_cast<std::size_t>(e.a)];
                    const int vb = tet.vertexIds[static_cast<std::size_t>(e.b)];
                    const auto eit = edgeIndexByKey.find(edgeKey(va, vb));
                    if (eit != edgeIndexByKey.end()) {
                        // Determine the within-edge sub-index k in [0, order-1).
                        // table.lambda[n] has lambda_a = 1 - t and lambda_b = t
                        // with t = (k+1)/order; we want k.
                        const double t = table.lambda[n][static_cast<std::size_t>(e.b)];
                        int subIdx = static_cast<int>(std::lround(t * fieldOutputOrder)) - 1;
                        if (subIdx < 0) subIdx = 0;
                        if (subIdx >= edgePointsPerEdge) subIdx = edgePointsPerEdge - 1;
                        // Order the two interior points on each edge by the
                        // *lower-vertex-id-first* direction, so neighbouring
                        // cells that traverse the edge in opposite orientation
                        // still hit the same global id.
                        if (va > vb) {
                            subIdx = (edgePointsPerEdge - 1) - subIdx;
                        }
                        globalId = edgeBase + eit->second * edgePointsPerEdge + subIdx;
                    }
                } else if (table.localFaceCenter[n] >= 0 && faceInteriorPerFace > 0) {
                    const auto& f = kVtkFaces[static_cast<std::size_t>(table.localFaceCenter[n])];
                    const std::array<int, 3> tri = faceKey(
                        tet.vertexIds[static_cast<std::size_t>(f.a)],
                        tet.vertexIds[static_cast<std::size_t>(f.b)],
                        tet.vertexIds[static_cast<std::size_t>(f.c)]);
                    const auto fit = faceSlotByKey.find(tri);
                    if (fit != faceSlotByKey.end()) {
                        globalId = faceBase + fit->second;
                    }
                } else {
                    // Volume interior: unique to this cell.
                    globalId = volumeCursor++;
                }
                cellLocalToGlobal[cellIndex][n] = globalId;
            }
        }
    }

    // Reconstruct per-cell field samples (cell-local), then average into
    // global slots based on cellLocalToGlobal.
    const auto eFieldLocal = reconstructLagrangePointField(mesh, topology, edgeDofs, table);
    const auto hFieldLocal = reconstructLagrangeMagneticField(mesh, project, topology, edgeDofs, table, frequencyHz);

    std::vector<Vec3> globalPoints(totalGlobalNodes, Vec3{0.0, 0.0, 0.0});
    std::vector<ComplexVec3> eField(totalGlobalNodes, ComplexVec3{});
    std::vector<ComplexVec3> hField(totalGlobalNodes, ComplexVec3{});
    std::vector<int> incidence(totalGlobalNodes, 0);

    for (std::size_t cellIndex = 0; cellIndex < totalCells; ++cellIndex) {
        const auto& tet = mesh.tetrahedra[cellIndex];
        std::array<Vec3, 4> p{};
        bool valid = true;
        for (int i = 0; i < 4; ++i) {
            const auto pit = mesh.pointsById.find(tet.vertexIds[static_cast<std::size_t>(i)]);
            if (pit == mesh.pointsById.end()) {
                valid = false;
                break;
            }
            p[static_cast<std::size_t>(i)] = pit->second;
        }
        if (!valid) {
            continue;
        }
        const std::size_t base = cellIndex * nodesPerCell;
        for (std::size_t n = 0; n < nodesPerCell; ++n) {
            const int g = cellLocalToGlobal[cellIndex][n];
            if (g < 0) {
                continue;
            }
            const auto& lambda = table.lambda[n];
            const Vec3 xyz = lambda[0] * p[0] + lambda[1] * p[1]
                           + lambda[2] * p[2] + lambda[3] * p[3];
            // Geometry: same xyz from any incident cell, take any.
            globalPoints[static_cast<std::size_t>(g)] = xyz;
            // Field: accumulate, average at end.
            for (int c = 0; c < 3; ++c) {
                eField[static_cast<std::size_t>(g)][c] += eFieldLocal[base + n][c];
                hField[static_cast<std::size_t>(g)][c] += hFieldLocal[base + n][c];
            }
            ++incidence[static_cast<std::size_t>(g)];
        }
    }

    // Average shared values. Volume-interior nodes have incidence = 1 (no-op
    // divide); shared nodes get (sum of values) / (cell count). Tangential
    // continuity makes the sum already coherent in the tangential plane;
    // the divide handles the small normal mismatch as an in-place average.
    for (int g = 0; g < totalGlobalNodes; ++g) {
        const int c = incidence[static_cast<std::size_t>(g)];
        if (c <= 1) continue;
        const double inv = 1.0 / static_cast<double>(c);
        for (int k = 0; k < 3; ++k) {
            eField[static_cast<std::size_t>(g)][k] *= inv;
            hField[static_cast<std::size_t>(g)][k] *= inv;
        }
    }

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"2.2\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << totalGlobalNodes << "\" NumberOfCells=\"" << totalCells << "\">\n";
    out << "      <PointData Vectors=\"E_real\">\n";
    writeVectorDataArray(out, "E_real", eField, [](const auto& z) { return z.real(); });
    writeVectorDataArray(out, "E_imag", eField, [](const auto& z) { return z.imag(); });
    writeVectorDataArray(out, "H_real", hField, [](const auto& z) { return z.real(); });
    writeVectorDataArray(out, "H_imag", hField, [](const auto& z) { return z.imag(); });
    out << "      </PointData>\n";
    out << "      <CellData>\n";
    out << "      </CellData>\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    out << std::setprecision(15);
    for (const auto& xyz : globalPoints) {
        out << "          " << xyz.x << ' ' << xyz.y << ' ' << xyz.z << '\n';
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";
    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t cellIndex = 0; cellIndex < totalCells; ++cellIndex) {
        out << "         ";
        for (std::size_t n = 0; n < nodesPerCell; ++n) {
            out << ' ' << cellLocalToGlobal[cellIndex][n];
        }
        out << '\n';
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 1; i <= totalCells; ++i) {
        out << "          " << i * nodesPerCell << '\n';
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
#ifdef BPFEM_OUTPUT_PARAVIEW
    // ParaView path: VTK Lagrange Tetrahedron (cell type 71) carrying
    // HigherOrderDegrees. Renders correctly in ParaView >= 5.5 / VTK >= 8.2;
    // VisIt does not understand this cell type.
    //   10 = VTK_TETRA               (linear, 4 nodes)
    //   71 = VTK_LAGRANGE_TETRAHEDRON (arbitrary order; carries HigherOrderDegrees)
    const unsigned int cellType = (fieldOutputOrder == 1) ? 10U : 71U;
    for (std::size_t i = 0; i < totalCells; ++i) {
        out << "          " << cellType << '\n';
    }
    out << "        </DataArray>\n";
    if (fieldOutputOrder > 1) {
        out << "        <DataArray type=\"Int32\" Name=\"HigherOrderDegrees\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (std::size_t i = 0; i < totalCells; ++i) {
            out << "          " << fieldOutputOrder << ' ' << fieldOutputOrder << ' ' << fieldOutputOrder << '\n';
        }
        out << "        </DataArray>\n";
    }
#else
    // Default (VisIt-friendly) path:
    //   order 1 -> 10 = VTK_TETRA            (linear, 4 nodes)
    //   order 2 -> 24 = VTK_QUADRATIC_TETRA  (4 corners + 6 edge midpoints,
    //                                         in the same edge order our
    //                                         Lagrange table emits)
    // No HigherOrderDegrees array is written: cell type 24 has a fixed layout.
    const unsigned int cellType = (fieldOutputOrder == 1) ? 10U : 24U;
    for (std::size_t i = 0; i < totalCells; ++i) {
        out << "          " << cellType << '\n';
    }
    out << "        </DataArray>\n";
#endif
    out << "      </Cells>\n";
    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
}

std::vector<OutputWriter::ComplexVec3> OutputWriter::reconstructLagrangePointField(
    const Mesh& mesh,
    const EdgeTopology& topology,
    const std::vector<std::complex<double>>& edgeDofs,
    const LagrangeTetTable& table) {
    const std::size_t nodesPerCell = table.lambda.size();
    std::vector<ComplexVec3> values(mesh.tetrahedra.size() * nodesPerCell);
    const auto& elementDofs = topology.elementDofs();

    for (std::size_t elementIndex = 0; elementIndex < mesh.tetrahedra.size(); ++elementIndex) {
        const auto& tet = mesh.tetrahedra[elementIndex];
        std::array<Vec3, 4> p{};
        bool valid = true;
        for (int i = 0; i < 4; ++i) {
            const auto pit = mesh.pointsById.find(tet.vertexIds[static_cast<std::size_t>(i)]);
            if (pit == mesh.pointsById.end()) {
                valid = false;
                break;
            }
            p[static_cast<std::size_t>(i)] = pit->second;
        }
        if (!valid || tetraVolume(p[0], p[1], p[2], p[3]) <= 1.0e-24) {
            continue;
        }

        const auto grad = tetraGradients(p);
        const std::size_t base = elementIndex * nodesPerCell;
        for (std::size_t n = 0; n < nodesPerCell; ++n) {
            const auto& lambda = table.lambda[n];
            ComplexVec3 cell{};
            for (const auto& ref : elementDofs[elementIndex]) {
                if (ref.globalIndex < 0 || static_cast<std::size_t>(ref.globalIndex) >= edgeDofs.size()) {
                    continue;
                }
                const BasisValue basis = evaluateBasis(ref, lambda, grad);
                const std::complex<double> coeff = edgeDofs[static_cast<std::size_t>(ref.globalIndex)];
                cell[0] += coeff * basis.value.x;
                cell[1] += coeff * basis.value.y;
                cell[2] += coeff * basis.value.z;
            }
            values[base + n] = cell;
        }
    }

    return values;
}

std::vector<OutputWriter::ComplexVec3> OutputWriter::reconstructLagrangeMagneticField(
    const Mesh& mesh,
    const ProjectDefinition& project,
    const EdgeTopology& topology,
    const std::vector<std::complex<double>>& edgeDofs,
    const LagrangeTetTable& table,
    double frequencyHz) {
    const std::size_t nodesPerCell = table.lambda.size();
    std::vector<ComplexVec3> values(mesh.tetrahedra.size() * nodesPerCell);
    const auto& elementDofs = topology.elementDofs();
    const double omega = 2.0 * pi * frequencyHz;
    if (omega <= 0.0) {
        return values;
    }
    const std::complex<double> jUnit(0.0, 1.0);

    for (std::size_t elementIndex = 0; elementIndex < mesh.tetrahedra.size(); ++elementIndex) {
        const auto& tet = mesh.tetrahedra[elementIndex];
        std::array<Vec3, 4> p{};
        bool valid = true;
        for (int i = 0; i < 4; ++i) {
            const auto pit = mesh.pointsById.find(tet.vertexIds[static_cast<std::size_t>(i)]);
            if (pit == mesh.pointsById.end()) {
                valid = false;
                break;
            }
            p[static_cast<std::size_t>(i)] = pit->second;
        }
        if (!valid || tetraVolume(p[0], p[1], p[2], p[3]) <= 1.0e-24) {
            continue;
        }

        const auto grad = tetraGradients(p);
        const Material material = materialForBody(mesh, project, tet.bodyId);
        const double muR = std::max(material.relativePermeability, 1.0e-30);
        // Faraday law (e^{+j omega t} convention): -j omega mu_r mu_0 H = curl(E)
        // => H = j curl(E) / (omega mu_r mu_0)
        const std::complex<double> hScale = jUnit / (omega * muR * mu0);
        const std::size_t base = elementIndex * nodesPerCell;
        for (std::size_t n = 0; n < nodesPerCell; ++n) {
            const auto& lambda = table.lambda[n];
            ComplexVec3 curlE{};
            for (const auto& ref : elementDofs[elementIndex]) {
                if (ref.globalIndex < 0 || static_cast<std::size_t>(ref.globalIndex) >= edgeDofs.size()) {
                    continue;
                }
                const BasisValue basis = evaluateBasis(ref, lambda, grad);
                const std::complex<double> coeff = edgeDofs[static_cast<std::size_t>(ref.globalIndex)];
                curlE[0] += coeff * basis.curl.x;
                curlE[1] += coeff * basis.curl.y;
                curlE[2] += coeff * basis.curl.z;
            }
            values[base + n][0] = hScale * curlE[0];
            values[base + n][1] = hScale * curlE[1];
            values[base + n][2] = hScale * curlE[2];
        }
    }

    return values;
}

}  // namespace fem
