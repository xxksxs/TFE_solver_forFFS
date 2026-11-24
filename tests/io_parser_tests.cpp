#include "bpfem/io/AEDTParser.hpp"
#include "bpfem/io/NGMeshParser.hpp"
#include "bpfem/io/PortFaceResolver.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/bc/WavePortBC.hpp"
#include "bpfem/core/Constants.hpp"
#include "bpfem/core/Utilities.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"
#include "bpfem/tfe/TransfiniteElementBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string_view>

namespace {

// 条件不满足时让测试立即失败并给出简短原因。
void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

// 在临时目录写入最小 AEDT/NGMesh 文本，用来覆盖 Objects(...) 端口路径。
void writeFixtures(const std::filesystem::path& aedtPath,
                   const std::filesystem::path& meshPath) {
    {
        std::ofstream out(aedtPath);
        out << "$begin '1'\n"
            << "ID=0\n"
            << "BoundType='Wave Port'\n"
            << "Objects(42)\n"
            << "NumModes=1\n"
            << "$end '1'\n"
            << "Frequency='95GHz'\n"
            << "RangeStart='90GHz'\n"
            << "RangeEnd='100GHz'\n"
            << "RangeCount=101\n"
            << "SourceEntry(ID=0, Index=0, Terminal=false, Terminated=false, "
               "Magnitude='1W', Phase='0deg')\n";
    }
    {
        std::ofstream out(meshPath);
        out << "user_unit_name mm\n"
            << "user_units_per_one_meter 1000\n"
            << "body_id 42 body_name PortSheet nvelems_on_body 0\n"
            << "face_ids 99\n"
            << "bbox_xmin 0\n"
            << "bbox_ymin 0\n"
            << "bbox_zmin 0\n"
            << "bbox_xmax 1\n"
            << "bbox_ymax 1\n"
            << "bbox_zmax 0\n"
            << "pid 1 coords 0 0 0\n"
            << "pid 2 coords 1 0 0\n"
            << "pid 3 coords 0 1 0\n"
            << "pid 4 coords 0 0 1\n"
            << "facet_id 5 face_ids 7\n"
            << "seid 1 facet_id 5 vert_ids 1 2 3\n"
            << "veid 1 body_id 6 vert_ids 1 2 3 4\n";
    }
}

// 验证 AEDT objectId、NGMesh body 包围盒和几何面解析能够串联工作。
void testObjectPortResolution() {
    const auto tempDir = std::filesystem::temp_directory_path() / "bpfem_io_parser_tests";
    std::filesystem::create_directories(tempDir);
    const auto aedtPath = tempDir / "fixture.aedt";
    const auto meshPath = tempDir / "fixture.ngmesh";
    writeFixtures(aedtPath, meshPath);

    fem::ProjectDefinition project = fem::AEDTParser{}.parse(aedtPath);
    const fem::Mesh mesh = fem::NGMeshParser{}.parse(meshPath);
    require(project.ports.size() == 1, "AEDT object port count mismatch");
    require(project.ports[0].objectId == 42, "AEDT objectId mismatch");
    require(project.ports[0].faceId == -1, "Object port should be unresolved before mesh matching");
    require(project.ports[0].excited, "AEDT source excitation mismatch");
    require(mesh.bodies.at(42).hasBounds, "NGMesh body bounds were not parsed");
    require(mesh.bodies.at(42).faceIds.size() == 1, "NGMesh standalone face_ids were not parsed");

    fem::PortFaceResolver::resolve(project, mesh);
    require(project.ports[0].faceId == 7, "Object port resolved to the wrong mesh face");

    std::filesystem::remove_all(tempDir);
}

void testDielectricInterfaceAssembly() {
    const auto path = std::filesystem::temp_directory_path() / "bpfem_dielectric_interface.aedt";
    {
        std::ofstream out(path);
        out << "$begin 'Materials'\n$begin 'ceramic'\npermittivity='4'\n"
            << "thermal_conductivity='25'\n$end 'ceramic'\n$end 'Materials'\n"
            << "$begin 'Attributes'\nName='Window_1'\nMaterialValue='\"ceramic\"'\n$end 'Attributes'\n"
            << "$begin '5'\nBoundType='Perfect E'\nFaces(10, 13)\n$end '5'\n";
    }
    auto project = fem::AEDTParser{}.parse(path);
    std::filesystem::remove(path);
    require(project.bodyMaterials.at("Window_1") == "ceramic", "AEDT object material mapping lost");
    require(project.materials.at("ceramic").conductivity == 0.0,
            "Thermal conductivity was used as electrical conductivity");
    require(project.explicitPecFaces == std::vector<int>({10, 13}), "AEDT Perfect E faces were lost");
    auto pec = project;
    project.explicitPecFaces.clear();

    fem::Mesh mesh;
    mesh.pointsById = {{1, {0.0, 0.0, 0.0}}, {2, {1.0, 0.0, 0.0}}, {3, {0.0, 1.0, 0.0}},
                       {4, {0.0, 0.0, 1.0}}, {5, {0.0, 0.0, -1.0}}};
    mesh.bodies[6].name = "vacuum";
    mesh.bodies[40].name = "Window_1";
    mesh.tetrahedra = {{1, 6, {1, 2, 3, 4}}, {2, 40, {1, 2, 3, 5}}};
    mesh.facets[1].faceId = 10;
    mesh.facets[2].faceId = 11;
    mesh.surfaceTriangles = {{1, 1, {1, 2, 3}}, {2, 2, {1, 2, 4}}};
    for (int order : {0, 1}) {
        fem::EdgeTopology topology(mesh, order);
        fem::PortModeSolver modes(mesh, topology);
        fem::FEMAssembler mixed(mesh, project, topology, modes);
        const auto affine = mixed.buildAffineSystem();
        const int interior = topology.edgeZeroDofForNodes(1, 3);
        const int exterior = topology.edgeZeroDofForNodes(1, 2);
        require(!affine.constrainedDofs.count(interior), "Material interface was treated as PEC");
        require(affine.constrainedDofs.count(exterior) == 1, "Exterior PEC edge became unconstrained");
        if (order == 1) {
            for (const int dof : topology.faceFirstDofsForNodes(1, 2, 3)) {
                require(!affine.constrainedDofs.count(dof), "High-order interface trace was treated as PEC");
            }
        }
        auto homogeneous = project;
        homogeneous.bodyMaterials.clear();
        const auto vacuum = fem::FEMAssembler(mesh, homogeneous, topology, modes).buildAffineSystem();
        std::vector<std::complex<double>> x(affine.M.size());
        x[interior] = 1.0;
        const auto mixedMass = affine.M.multiply(x)[interior];
        const auto vacuumMass = vacuum.M.multiply(x)[interior];
        require(std::abs(mixedMass - 2.5 * vacuumMass) <= 1e-12 * std::abs(vacuumMass),
                "Dielectric was not applied to its body mass contribution");
        const auto internalWall = fem::FEMAssembler(mesh, pec, topology, modes).buildAffineSystem();
        require(internalWall.constrainedDofs.count(interior) == 1, "Explicit internal PEC lost its constraint");
    }
}

// 对用户指定的真实工程执行纯解析检查，不建立拓扑也不启动 FEM 求解。
void inspectRealInputs(const std::filesystem::path& aedtPath,
                       const std::filesystem::path& meshPath) {
    fem::ProjectDefinition project = fem::AEDTParser{}.parse(aedtPath);
    const fem::Mesh mesh = fem::NGMeshParser{}.parse(meshPath);
    fem::PortFaceResolver::resolve(project, mesh);

    std::cout << "points=" << mesh.pointsById.size()
              << " surface_triangles=" << mesh.surfaceTriangles.size()
              << " tetrahedra=" << mesh.tetrahedra.size() << '\n';
    for (const auto& port : project.ports) {
        std::cout << "port=" << port.id
                  << " object=" << port.objectId
                  << " face=" << port.faceId
                  << " modes=" << port.modes
                  << " excited=" << (port.excited ? "true" : "false") << '\n';
    }
}

// Run on the same parsed inputs and frequencies as the Direct smoke run.
// Keep this opt-in so the small no-argument CTest needs no external cases.
void inspectPortBaseline(const std::filesystem::path& aedtPath,
                         const std::filesystem::path& meshPath, int modesPerPort) {
    using Builder = fem::fastsweep::PolynomialPortMomentBuilder;
    using Complex = Builder::Complex;
    using Vector = std::vector<Complex>;
    using fem::fastsweep::bilinear;
    using fem::fastsweep::norm2;
    constexpr double tolerance = 1e-11;
    auto project = fem::AEDTParser{}.parse(aedtPath);
    const auto mesh = fem::NGMeshParser{}.parse(meshPath);
    fem::PortFaceResolver::resolve(project, mesh);
    require(project.ports.size() == 2, "The port baseline requires two physical ports");
    fem::EdgeTopology topology(mesh, 0);
    fem::PortModeSolver modes(mesh, topology);
    if (modesPerPort == 2) {
        fem::tfe::TfeOptions options;
        options.modesPerPort = modesPerPort;
        const fem::tfe::TransfiniteElementBuilder builder(mesh, topology, modes, options);
        for (const auto& port : project.ports) modes.setMultiMode(port.faceId, builder.build(port.faceId));
    }

    fem::FEMAssembler assembler(mesh, project, topology, modes);
    const auto wavePort = std::make_shared<fem::bc::WavePortBC>();
    assembler.setBoundaryConditions({wavePort});
    const auto affine = assembler.buildAffineSystem();
    const auto ports = Builder::buildPortVectors(affine, affine.K.size());
    require(affine.lossless, "The baseline requires a lossless assembled device");
    require(ports.size() == project.ports.size() * static_cast<std::size_t>(modesPerPort),
            "Wrong number of virtual ports");
    require(std::count(affine.isExcitationMode.begin(), affine.isExcitationMode.end(), true) == 1,
            "Wrong number of input modes");
    require(affine.isExcitationMode.front(), "Expected to excite only the first physical port, first mode");
    const auto frequencies = fem::buildFrequencies(project.sweep, project.sweep.count);
    require(!frequencies.empty(), "Empty baseline frequency list");

    Vector x(affine.K.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i] = {std::sin(0.77 * i), std::cos(0.49 * i)};
    }
    std::cout << std::setprecision(17) << "{\"matrix_tolerance\":" << tolerance
              << ",\"pec_dofs\":" << affine.constrainedDofs.size() << ",\"channels\":[";
    for (std::size_t v = 0; v < ports.size(); ++v) {
        if (v != 0) std::cout << ',';
        const auto portIndex = static_cast<std::size_t>(affine.projectPortIndex.at(v));
        const auto modeIndex = v % static_cast<std::size_t>(modesPerPort);
        const auto& port = project.ports.at(portIndex);
        const auto* multi = modes.multiMode(port.faceId);
        const auto& mode = multi ? multi->modes.at(modeIndex) : modes.solve(port.faceId);
        std::vector<std::pair<int, double>> trace;
        for (const auto& entry : mode.couplingWeights) {
            if (!affine.constrainedDofs.count(entry.first)) trace.push_back(entry);
        }
        require(!trace.empty() && trace == affine.portCoupling.at(v), "Missing or misrouted modal trace");
        require(mode.cutoffWavenumberSquared == affine.portCutoffSquared.at(v), "Wrong modal cutoff mapping");
        const double cutoffHz = fem::c0 * std::sqrt(mode.cutoffWavenumberSquared) / (2.0 * fem::pi);
        std::cout << "{\"port_id\":" << port.id << ",\"mode\":" << modeIndex + 1
                  << ",\"cutoff_hz\":" << cutoffHz << ",\"trace_entries\":" << trace.size()
                  << ",\"excited\":" << (affine.isExcitationMode[v] ? "true" : "false") << '}';
    }

    const auto relativeDifference = [](Vector a, const Vector& b) {
        for (std::size_t i = 0; i < a.size(); ++i) a[i] -= b[i];
        const double expectedNorm = norm2(b);
        return expectedNorm > 0.0 ? norm2(a) / expectedNorm : norm2(a);
    };
    std::cout << "],\"samples\":[";
    bool first = true;
    for (const auto f : frequencies) {
        const auto series = Builder::buildTaylorSeries(project, modes, affine, f, 3,
                                                       Builder::Variable::WavenumberShift);
        Vector directRhs;
        const auto direct = assembler.assemble(f, directRhs);
        const auto a0 = Builder::applyMatrixCoefficient(affine, series, ports, 0, x);
        const auto b0 = Builder::rhsCoefficient(series, ports, x.size(), 0);
        const double matrixError = relativeDifference(a0, direct.multiply(x));
        const double rhsError = relativeDifference(b0, directRhs);
        require(matrixError <= tolerance && rhsError <= tolerance && norm2(directRhs) > 0.0,
                "Direct and Taylor zero-order systems differ");

        // Compare the actual boundary block separately, so stiffness scaling
        // cannot hide a bad mode contribution in the full-system residual.
        fem::SparseMatrix boundary(x.size());
        Vector boundaryRhs(x.size()), expectedBoundary(x.size());
        wavePort->apply(boundary, boundaryRhs, {f, mesh, project, topology, modes, affine.constrainedDofs});
        for (std::size_t v = 0; v < ports.size(); ++v) {
            const auto factor = Complex(0.0, series.admittance[v][0]) * bilinear(ports[v], x);
            for (const auto& [dof, weight] : affine.portCoupling[v]) expectedBoundary[dof] += factor * weight;
            if (!affine.isExcitationMode[v]) {
                require(series.excitation[v][0] == Complex{}, "Passive mode received an incident RHS");
            }
        }
        const double boundaryError = relativeDifference(expectedBoundary, boundary.multiply(x));
        const double boundaryRhsError = relativeDifference(b0, boundaryRhs);
        require(boundaryError <= tolerance && boundaryRhsError <= tolerance, "Wrong assembled wave-port block");
        const auto a1 = Builder::applyMatrixCoefficient(affine, series, ports, 1, x);
        const auto b1 = Builder::rhsCoefficient(series, ports, x.size(), 1);
        for (const int dof : affine.constrainedDofs) {
            require(a0[dof] == x[dof] && a1[dof] == Complex{} && b0[dof] == Complex{} && b1[dof] == Complex{},
                    "PEC row leaked into Taylor derivatives or excitation");
        }
        if (!first) std::cout << ',';
        first = false;
        std::cout << "{\"frequency_hz\":" << f << ",\"a0_relative_error\":" << matrixError
                  << ",\"b0_relative_error\":" << rhsError << ",\"boundary_relative_error\":" << boundaryError
                  << ",\"boundary_rhs_relative_error\":" << boundaryRhsError << '}';
    }
    std::cout << "]}\n";
}

}  // namespace

// 无参数时运行 CTest；给定 AEDT 和 NGMesh 时执行真实输入的轻量映射检查。
int main(int argc, char** argv) {
    if (argc == 5 && std::string_view(argv[1]) == "--port-baseline") {
        try {
            const std::string_view mode(argv[4]);
            require(mode == "numerical" || mode == "tfe", "Expected numerical or tfe");
            inspectPortBaseline(argv[2], argv[3], mode == "tfe" ? 2 : 1);
            return 0;
        } catch (const std::exception& error) {
            std::cerr << error.what() << '\n';
            return 1;
        }
    }
    if (argc == 3) {
        inspectRealInputs(argv[1], argv[2]);
        return 0;
    }
    testObjectPortResolution();
    testDielectricInterfaceAssembly();
    return 0;
}
