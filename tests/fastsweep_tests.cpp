#include "bpfem/fastsweep/GalerkinReducedModel.hpp"
#include "bpfem/fastsweep/PadeApproximant.hpp"
#include "bpfem/fastsweep/PolynomialMomentRecurrence.hpp"
#include "bpfem/fastsweep/WellConditionedBasisBuilder.hpp"
#include "bpfem/fem/EdgeTopology.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Complex = std::complex<double>;

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void requireNear(double actual, double expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message);
    }
}

void requireNear(Complex actual, Complex expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message);
    }
}

void testPolynomialMomentRecurrence() {
    std::vector<std::vector<Complex>> rhs = {
        {Complex(1.0, 0.0)},
        {Complex(0.0, 0.0)},
        {Complex(0.0, 0.0)}
    };
    std::vector<fem::fastsweep::PolynomialMomentRecurrence::LinearOperator> ops;
    ops.push_back([](const std::vector<Complex>& x) {
        return std::vector<Complex>{x[0]};
    });
    auto solve = [](const std::vector<Complex>& b) {
        return std::vector<Complex>{b[0] / Complex(2.0, 0.0)};
    };

    const auto moments =
        fem::fastsweep::PolynomialMomentRecurrence::generatePolynomialMoments(rhs, ops, solve);
    require(moments.size() == 3, "PolynomialMomentRecurrence returned wrong count");
    requireNear(moments[0][0], Complex(0.5, 0.0), 1.0e-14, "moment 0 mismatch");
    requireNear(moments[1][0], Complex(-0.25, 0.0), 1.0e-14, "moment 1 mismatch");
    requireNear(moments[2][0], Complex(0.125, 0.0), 1.0e-14, "moment 2 mismatch");
}

void testPadeApproximant() {
    const std::vector<Complex> moments = {
        Complex(1.0, 0.0),
        Complex(2.0, 0.0),
        Complex(4.0, 0.0)
    };
    const auto pade = fem::fastsweep::PadeApproximant::build(moments, 0, 1);
    require(pade.numeratorOrder() == 0, "Pade numerator order mismatch");
    require(pade.denominatorOrder() == 1, "Pade denominator order mismatch");
    requireNear(pade.denominator()[1], Complex(-2.0, 0.0), 1.0e-14,
                "Pade denominator coefficient mismatch");
    requireNear(pade.evaluate(Complex(0.1, 0.0)), Complex(1.25, 0.0), 1.0e-14,
                "Pade evaluation mismatch");
}

void testWellConditionedBasisBuilder() {
    fem::fastsweep::WellConditionedBasisBuilder builder(1.0e-14);
    require(builder.append({Complex(1.0, 0.0), Complex(0.0, 0.0)}),
            "WCAWE builder rejected first vector");
    require(builder.append({Complex(1.0, 0.0), Complex(1.0e-8, 0.0)}),
            "WCAWE builder rejected independent near-collinear vector");
    require(builder.dimension() == 2, "WCAWE builder dimension mismatch");
    require(builder.conditionRecords().size() == 2, "WCAWE condition record count mismatch");
    require(builder.conditionRecords().back().aweConditionProxy > 1.0e6,
            "WCAWE condition proxy did not detect ill-conditioned moment basis");
    require(builder.conditionRecords().back().orthogonalityError < 1.0e-10,
            "WCAWE orthogonality error too large");
    require(builder.maxMomentReconstructionError() < 1.0e-10,
            "WCAWE moment reconstruction error too large");
}

fem::PortMode makePortMode(int faceId, int dof) {
    fem::PortMode mode;
    mode.faceId = faceId;
    mode.cutoffWavenumberSquared = 0.0;
    mode.couplingWeights.push_back({dof, 1.0});
    fem::PortQuadraturePoint qp;
    qp.modeFieldValue = {1.0, 0.0, 0.0};
    qp.normal = {0.0, 0.0, 1.0};
    qp.weight = 1.0;
    mode.quadrature.push_back(qp);
    return mode;
}

void testGalerkinReducedModel() {
    fem::ProjectDefinition project;
    project.ports.push_back({1, 10, 1, true, 1.0, 0.0});
    project.ports.push_back({2, 20, 1, false, 0.0, 0.0});

    fem::FEMAssembler::AffineSystem affine;
    affine.K = fem::SparseMatrix(2);
    affine.K.add(0, 0, Complex(1.0, 0.0));
    affine.K.add(1, 1, Complex(1.0, 0.0));
    affine.M = fem::SparseMatrix(2);
    affine.portCoupling = {{{0, 1.0}}, {{1, 1.0}}};
    affine.portCutoffSquared = {0.0, 0.0};
    affine.portFaceIds = {10, 20};
    affine.projectPortIndex = {0, 1};
    affine.isExcitationMode = {true, false};
    affine.lossless = true;

    fem::Mesh mesh;
    fem::EdgeTopology topology(mesh, 0);
    fem::PortModeSolver portModeSolver(mesh, topology);
    portModeSolver.setPrecomputed(10, makePortMode(10, 0));
    portModeSolver.setPrecomputed(20, makePortMode(20, 1));

    std::vector<std::vector<Complex>> basis = {
        {Complex(1.0, 0.0), Complex(0.0, 0.0)},
        {Complex(0.0, 0.0), Complex(1.0, 0.0)}
    };
    std::vector<std::vector<Complex>> portVectors = {
        {Complex(1.0, 0.0), Complex(0.0, 0.0)},
        {Complex(0.0, 0.0), Complex(1.0, 0.0)}
    };

    fem::fastsweep::GalerkinReducedModel model;
    model.build(project, affine, portModeSolver, std::move(basis), portVectors);
    require(model.dimension() == 2, "Galerkin ROM dimension mismatch");
    require(model.basisOrthogonalityError() < 1.0e-14,
            "Galerkin basis orthogonality error too large");
    const auto sp = model.evaluate(1.0e9);
    require(std::isfinite(sp.s11.real()) && std::isfinite(sp.s11.imag()),
            "Galerkin S11 is not finite");
    require(std::isfinite(sp.s21.real()) && std::isfinite(sp.s21.imag()),
            "Galerkin S21 is not finite");
    const auto field = model.reconstructField(1.0e9);
    require(field.size() == 2, "Galerkin reconstructed field size mismatch");
}

}  // namespace

int main() {
    testPolynomialMomentRecurrence();
    testPadeApproximant();
    testWellConditionedBasisBuilder();
    testGalerkinReducedModel();
    return 0;
}
