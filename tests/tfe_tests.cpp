#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PolynomialMomentRecurrence.hpp"
#include "bpfem/fastsweep/TfePolynomialModel.hpp"
#include "bpfem/linalg/FactorizedSolveSession.hpp"
#include "bpfem/linalg/PardisoBackend.hpp"
#include "bpfem/tfe/TfeSystem.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
using namespace fem;
using fem::tfe::TfeSystem;
using Vector = TfeSystem::Vector;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

Mesh waveguide(int transverse, int longitudinal) {
    Mesh mesh;
    const auto node = [=](int x, int y, int z) {
        return 1 + x + (transverse + 1) * (y + (transverse + 1) * z);
    };
    for (int z = 0; z <= longitudinal; ++z) {
        for (int y = 0; y <= transverse; ++y) {
            for (int x = 0; x <= transverse; ++x) {
                mesh.pointsById.emplace(node(x, y, z),
                    Vec3{0.02 * x / transverse, 0.012 * y / transverse, 0.02 * z / longitudinal});
            }
        }
    }
    for (int z = 0; z < longitudinal; ++z) {
        for (int y = 0; y < transverse; ++y) {
            for (int x = 0; x < transverse; ++x) {
                std::array<int, 3> permutation{0, 1, 2};
                do {
                    std::array<int, 3> p{x, y, z};
                    Tetrahedron tet;
                    tet.id = static_cast<int>(mesh.tetrahedra.size());
                    tet.bodyId = 1;
                    tet.vertexIds[0] = node(x, y, z);
                    for (int i = 0; i < 3; ++i) {
                        ++p[permutation[i]];
                        tet.vertexIds[i + 1] = node(p[0], p[1], p[2]);
                    }
                    mesh.tetrahedra.push_back(tet);
                } while (std::next_permutation(permutation.begin(), permutation.end()));
            }
        }
    }
    std::map<std::array<int, 3>, int> faceCounts;
    for (const auto& tet : mesh.tetrahedra) {
        for (const auto& local : EdgeTopology::localFaces) {
            std::array<int, 3> face{tet.vertexIds[local[0]], tet.vertexIds[local[1]], tet.vertexIds[local[2]]};
            std::sort(face.begin(), face.end());
            ++faceCounts[face];
        }
    }
    for (const auto& [face, count] : faceCounts) {
        if (count != 1) continue;
        const bool entrance = std::all_of(face.begin(), face.end(),
            [&](int i) { return mesh.pointsById.at(i).z == 0.0; });
        const bool exit = std::all_of(face.begin(), face.end(),
            [&](int i) { return std::abs(mesh.pointsById.at(i).z - 0.02) < 1.0e-12; });
        const int id = entrance ? 1 : exit ? 2 : 3;
        mesh.facets[id].faceId = id;
        mesh.surfaceTriangles.push_back({static_cast<int>(mesh.surfaceTriangles.size()), id, face});
    }
    return mesh;
}

void validate(int basisOrder, int transverse, int longitudinal) {
    const Mesh mesh = waveguide(transverse, longitudinal);
    ProjectDefinition project;
    project.ports.resize(2);
    for (int i = 0; i < 2; ++i) {
        project.ports[i].id = i + 1;
        project.ports[i].faceId = i + 1;
        project.ports[i].modes = 6;
    }
    EdgeTopology topology(mesh, basisOrder);
    PortModeSolver modes(mesh, topology);
    FEMAssembler fem(mesh, project, topology, modes);
    TfeSystem system(mesh, project, topology, fem, modes, 0);
    require(system.channels().size() == 12, "Mode count mismatch");
    require(system.interiorCount() + system.traceCount() + fem.pecDofs().size() == topology.edgeCount(),
            "Invalid PEC/interior/trace partition");
    require(system.pattern()->n == system.interiorCount() + 12 && system.traceCount() > 12,
            "TFE did not replace trace DOFs");
    require(!system.channels().front().transverseMagnetic, "Fundamental mode must be TE");
    const auto& channels = system.channels();
    const auto tm = std::find_if(channels.begin(), channels.end(), [](const auto& c) {
        return c.portId == 1 && c.transverseMagnetic;
    });
    require(tm != channels.end(), "TM spectrum was omitted");
    const auto tmIndex = static_cast<std::size_t>(tm - channels.begin());
    const double te10 = c0 / 0.04;
    const double tm11 = c0 * std::hypot(1.0 / 0.02, 1.0 / 0.012) / 2.0;
    require(std::abs(channels.front().cutoffHz / te10 - 1.0) < 0.045, "TE10 cutoff error");
    require(std::abs(tm->cutoffHz / tm11 - 1.0) < 0.09, "TM11 cutoff error");
    require(system.admittance(0, te10 * 0.8).imag() < 0.0, "Wrong evanescent TE branch");
    require(system.admittance(tmIndex, tm11 * 0.8).imag() > 0.0, "Wrong evanescent TM branch");
    require(system.status(0, channels.front().cutoffHz) == tfe::ModeStatus::Cutoff, "Cutoff detection failed");

    constexpr double f = 9.5e9;
    const auto a = system.assemble(f);
    Vector incident(channels.size());
    incident[0] = 1.0;
    const auto b = system.rhs(f, incident);
    require(std::abs(b[system.interiorCount()] - std::complex<double>(0.0, 4.0 * pi * f * mu0)) < 1.0e-10,
            "Incorrect textbook incident factor");
    Vector z(a.size());
    for (std::size_t i = 0; i < z.size(); ++i) {
        z[i] = {std::sin(static_cast<double>(i)), std::cos(static_cast<double>(i))};
    }
    const auto edge = system.restoreEdgeField(f, z);
    for (const auto constrained : fem.pecDofs()) require(edge[constrained] == 0.0, "Nonzero PEC field");
    const auto vol = fem.buildVolumeSystem();
    const double k = 2.0 * pi * f / c0;
    const auto lhs = fastsweep::bilinear(z, a.multiply(z));
    auto expected = fastsweep::bilinear(edge, vol.K.multiply(edge))
                  - k * k * fastsweep::bilinear(edge, vol.M.multiply(edge));
    for (std::size_t i = system.interiorCount(); i < z.size(); ++i) {
        expected += std::complex<double>(0.0, k * c0 * mu0) * z[i] * z[i];
    }
    require(std::abs(lhs - expected) < 1.0e-10 * std::max(1.0, std::abs(expected)), "Incorrect Q^T D Q");

    fastsweep::TfePolynomialModel polynomial(system, 9.0e9, 10.0e9);
    for (const auto hz : polynomial.frequencies()) {
        const auto exact = system.assemble(hz).multiply(z);
        auto error = polynomial.apply(polynomial.weights(polynomial.coordinate(hz)), z);
        for (std::size_t i = 0; i < error.size(); ++i) error[i] -= exact[i];
        require(fastsweep::norm2(error) <= 1.0e-10 * fastsweep::norm2(exact), "Chebyshev node mismatch");
    }
    const double t = 0.25;
    const double offset = -0.4;
    const auto w = polynomial.weights(t);
    const auto d = polynomial.weights(t, 1);
    const auto dd = polynomial.weights(t, 2);
    const auto target = polynomial.weights(t + offset);
    for (std::size_t i = 0; i < w.size(); ++i) {
        require(std::abs(w[i] + offset * d[i] + offset * offset * dd[i] - target[i]) < 1.0e-14,
                "Quadratic coefficient shift failed");
    }
#ifdef BPFEM_USE_MKL
    linalg::PardisoBackend solver;
    const auto result = solver.solve(a, b);
    const auto outgoing = system.outgoing(result.field, incident);
    require(std::abs(outgoing[0]) < 0.08, "Uniform waveguide reflection too large");
    require(std::abs(std::norm(outgoing[0]) + std::norm(outgoing[6]) - 1.0) < 1.0e-8, "Power is not conserved");
    const double beta = std::sqrt(k * k - std::pow(pi / 0.02, 2));
    const auto reference = std::polar(1.0, -0.02 * beta);
    // The sign convention for the output eigenvector is arbitrary.
    require(std::min(std::abs(outgoing[6] - reference), std::abs(outgoing[6] + reference)) < 0.12,
            "Uniform waveguide transmission phase is incorrect");
    Vector reverseIncident(channels.size());
    reverseIncident[6] = 1.0;
    const auto reverse = solver.solve(a, system.rhs(f, reverseIncident));
    require(std::abs(system.outgoing(reverse.field, reverseIncident)[0] - outgoing[6]) < 1.0e-10,
            "Power-normalized modal transmission is not reciprocal");

    auto a0 = polynomial.matrix(polynomial.weights(0.0));
    linalg::FactorizedSolveSession session(solver, a0);
    const auto weights1 = polynomial.weights(0.0, 1);
    const auto weights2 = polynomial.weights(0.0, 2);
    std::vector<Vector> rhsCoefficients{b, b, b};
    for (auto& entry : rhsCoefficients[1]) entry *= 0.5e9 / f;
    // Artificial nonzero y2 checks the generic recurrence independently
    // of TFE's frequency-linear excitation.
    for (auto& entry : rhsCoefficients[2]) entry *= 0.01;
    const std::vector<fastsweep::PolynomialMomentRecurrence::LinearOperator> operators{
        [&](const auto& x) { return polynomial.apply(weights1, x); },
        [&](const auto& x) { return polynomial.apply(weights2, x); }};
    const auto moments = fastsweep::PolynomialMomentRecurrence::generatePolynomialMoments(
        rhsCoefficients, operators, [&](const auto& input) { return session.solve(input).field; });
    require(moments.size() == 3, "Unexpected moment count");
    auto residual = a0.multiply(moments[2]);
    const auto first = polynomial.apply(weights1, moments[1]);
    const auto second = polynomial.apply(weights2, moments[0]);
    for (std::size_t i = 0; i < residual.size(); ++i) {
        residual[i] += first[i] + second[i] - rhsCoefficients[2][i];
    }
    require(fastsweep::norm2(residual) < 1.0e-10 * fastsweep::norm2(b), "Second moment misses rhs coefficient");
#endif
}

void validateRadiation() {
    const Mesh mesh = waveguide(2, 2);
    ProjectDefinition project;
    project.ports.resize(2);
    project.ports[0].faceId = 1;
    project.ports[1].faceId = 2;
    project.radiationFaces = {3};
    EdgeTopology topology(mesh, 1);
    PortModeSolver modes(mesh, topology);
    FEMAssembler fem(mesh, project, topology, modes);
    require(!fem.pecDofs().empty(), "Legacy wall mask is missing");
    fem.enableRadiationFaces();
    require(fem.pecDofs().empty(), "ABC face remained PEC");
    const auto abc = fem.buildVolumeSystem().radiation;
    Vector z(abc.size(), 1.0);
    const auto form = fastsweep::hdot(z, abc.multiply(z));
    require(form.real() > 0.0 && form.imag() == 0.0, "ABC trace mass is not dissipative");
}
}  // namespace

int main() {
    try {
        validate(0, 5, 8);
        validate(1, 4, 6);
        validateRadiation();
        std::cout << "TFE modal spectrum, block assembly, and waveguide checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
