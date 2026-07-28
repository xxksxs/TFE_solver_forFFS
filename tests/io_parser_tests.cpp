#include "bpfem/io/AEDTParser.hpp"
#include "bpfem/io/NGMeshParser.hpp"
#include "bpfem/io/PortFaceResolver.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

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

}  // namespace

// 无参数时运行 CTest；给定 AEDT 和 NGMesh 时执行真实输入的轻量映射检查。
int main(int argc, char** argv) {
    if (argc == 3) {
        inspectRealInputs(argv[1], argv[2]);
        return 0;
    }
    testObjectPortResolution();
    return 0;
}
