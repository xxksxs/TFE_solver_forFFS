#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fem/EdgeTopology.hpp"

#include <array>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace fem {

class OutputWriter {
public:
    // VTK Lagrange tetrahedron output order. Solver basis order is independent
    // of visualization order; pick the latter to faithfully render the chosen
    // hierarchical basis:
    //   order=1 -> 4 nodes per tet (legacy linear VTK_TETRA-equivalent)
    //   order=2 -> 10 nodes (corners + edge midpoints), captures EdgeFirst bubbles
    //   order=3 -> 20 nodes (+ face midpoints), additionally captures FaceFirst bubbles
    static constexpr int kMinFieldOutputOrder = 1;
    static constexpr int kMaxFieldOutputOrder = 3;

    static void writeSParameters(const std::filesystem::path& path, const std::vector<SParameterPoint>& data);

    // Writes a VTU containing only the electric field (E_real, E_imag) and the
    // reconstructed magnetic field (H_real, H_imag). Magnitudes / dB / PEC
    // tangential diagnostics are no longer emitted; ParaView's Calculator filter
    // can derive any scalar from these four vector arrays. The frequency is
    // required to map curl(E) -> H via H = j curl(E) / (omega * mu_r * mu_0).
    static void writeVTU(const std::filesystem::path& path,
                         const Mesh& mesh,
                         const ProjectDefinition& project,
                         const EdgeTopology& topology,
                         const std::vector<std::complex<double>>& edgeDofs,
                         double frequencyHz,
                         int fieldOutputOrder = 1);

private:
    using ComplexVec3 = std::array<std::complex<double>, 3>;

    struct LagrangeTetTable {
        int order = 1;                                   // p
        std::vector<std::array<double, 4>> lambda;       // barycentric coords per Lagrange node
        std::vector<int> faceOnly;                       // -1 if interior/corner; otherwise local face index (0..3) where lambda_k==0
        std::vector<int> cornerIndex;                    // local corner index 0..3 if this Lagrange node is a corner; -1 otherwise
        // Sharing classification (populated by buildLagrangeTetTable):
        //   localEdge   in [0, 6) when this node lives on cell edge kVtkEdges[localEdge]
        //               (interior of that edge); -1 otherwise.
        //   localFaceCenter in [0, 4) when this node is a face-interior node of
        //               cell face kVtkFaces[localFaceCenter]; -1 otherwise.
        // Volume-interior nodes have all four classification fields = -1.
        // For order <= 2 only corners + edges can be set; order=3 also fills
        // localFaceCenter for the four face-interior points.
        std::vector<int> localEdge;
        std::vector<int> localFaceCenter;
    };

    static LagrangeTetTable buildLagrangeTetTable(int order);
    static std::vector<ComplexVec3> reconstructLagrangePointField(const Mesh& mesh,
                                                                  const EdgeTopology& topology,
                                                                  const std::vector<std::complex<double>>& edgeDofs,
                                                                  const LagrangeTetTable& table);
    static std::vector<ComplexVec3> reconstructLagrangeMagneticField(const Mesh& mesh,
                                                                     const ProjectDefinition& project,
                                                                     const EdgeTopology& topology,
                                                                     const std::vector<std::complex<double>>& edgeDofs,
                                                                     const LagrangeTetTable& table,
                                                                     double frequencyHz);

    template <typename Func>
    static void writeVectorDataArray(std::ofstream& out, const std::string& name, const std::vector<ComplexVec3>& field, Func func) {
        out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        out << std::setprecision(15);
        for (const auto& v : field) {
            out << "          " << func(v[0]) << ' ' << func(v[1]) << ' ' << func(v[2]) << '\n';
        }
        out << "        </DataArray>\n";
    }
};

}  // namespace fem
