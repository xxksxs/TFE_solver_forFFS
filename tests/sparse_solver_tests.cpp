#include "bpfem/linalg/BiCGStabBackend.hpp"
#include "bpfem/app/Application.hpp"
#include "bpfem/core/RunEnvironment.hpp"
#include "bpfem/core/RunReport.hpp"
#include "bpfem/linalg/FactorizedSolveSession.hpp"
#include "bpfem/linalg/GmresBackend.hpp"
#include "bpfem/linalg/IterativeOperator.hpp"
#include "bpfem/linalg/MumpsBackend.hpp"
#include "bpfem/linalg/PardisoBackend.hpp"
#include "bpfem/linalg/PreconSchwarz.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace fem;
using namespace fem::linalg;

namespace {

using Vector = std::vector<std::complex<double>>;

class CountingIdentity final : public IPreconditioner {
public:
    void setup(const SparseMatrix&) override { ++preparations; }
    void apply(const Vector& x, Vector& y) const override { y = x; }
    const char* name() const override { return "counting identity"; }
    int preparations = 0;
};

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void testToleranceReport() {
    Options options;
    options.tolerance = 1e-7;
    options.mumpsBlrTolerance = 1.23456789012345e-10;
    RunReport report;
    fillReportFromEnvironment(report, RunEnvironment{}, options);
    int checked = 0;
    for (const auto& [key, value] : report.options) {
        if (key == "tolerance") {
            require(std::stod(value) == options.tolerance, "Linear tolerance did not round trip");
            ++checked;
        } else if (key == "mumpsBlrTolerance") {
            require(std::stod(value) == options.mumpsBlrTolerance, "BLR tolerance did not round trip");
            ++checked;
        }
    }
    require(checked == 2, "Missing reported tolerance");
}

#ifdef BPFEM_USE_MKL
void testMixedPrecisionGmres(const SparseMatrix& matrix, const Vector& expected) {
    GmresBackend solver;
    solver.setPreconditioner(std::make_shared<PreconPardisoSingle>());
    const SolverConfig config{30, 1e-12};
    {
        FactorizedSolveSession prepared(solver, matrix, config);
        const auto rhs = matrix.multiply(expected);
        const auto results = prepared.solveBatch({rhs, rhs});
        for (const auto& result : results) {
            require(result.status == SolveStatus::Converged && result.residual <= config.tolerance,
                    "Mixed precision iteration did not meet the original matrix tolerance");
            require(result.iterations > 1, "Expected a correction to the single precision solve");
        }
    }
    auto next = matrix;
    next.add(0, 0, {0.5, 0.9});
    const auto result = solver.solve(next, next.multiply(expected), config);
    require(result.status == SolveStatus::Converged && result.residual <= config.tolerance,
            "Mixed precision preconditioner must rebuild for changed coefficients");
}
#endif

SparseMatrix testMatrix() {
    SparseMatrix matrix(4);
    matrix.add(0, 0, {5.0, 1.0});
    matrix.add(0, 1, {1.0, 2.0});
    matrix.add(1, 1, {-4.0, 0.3});
    matrix.add(1, 2, {0.7, -0.8});
    matrix.add(2, 2, {6.0, 0.4});
    matrix.add(2, 3, {-0.4, 0.6});
    matrix.add(3, 3, {-5.0, -0.5});
    return matrix;
}

const Vector expected{{1.0, 1.0}, {2.0, -1.0}, {-0.7, 0.0}, {0.2, 0.8}};
const SolverConfig strict{80, 1e-11};

void check(const SolveResult& result, const Vector& solution) {
    require(result.status == SolveStatus::Converged, "Solve did not converge");
    require(result.residual <= strict.tolerance, "Original-matrix residual is too large");
    require(result.field.size() == solution.size(), "Incorrect solution size");
    for (std::size_t i = 0; i < solution.size(); ++i) {
        require(std::abs(result.field[i] - solution[i]) < 1e-9, "Wrong complex solution");
    }
}

template<class Backend>
void testIteration(const SparseMatrix& matrix) {
    Backend solver;
    auto precon = std::make_shared<CountingIdentity>();
    solver.setPreconditioner(precon);
    const auto rhs = matrix.multiply(expected);
    const Vector zero(expected.size());
    {
        FactorizedSolveSession session(solver, matrix, strict);
        const auto results = session.solveBatch({rhs, zero, rhs});
        check(results.at(0), expected);
        check(results.at(1), zero);
        check(results.at(2), expected);
        require(results[1].iterations == 0, "Zero RHS should require no iterations");
        require(precon->preparations == 1, "Prepared solve repeated preconditioner setup");
    }
    auto next = matrix;
    next.add(0, 0, {0.5, 0.9});
    check(solver.solve(next, next.multiply(expected), strict), expected);
    require(precon->preparations == 2, "New matrix did not rebuild the preconditioner");
    const auto failure = solver.solve(matrix, rhs, {1, 1e-15});
    require(failure.status == SolveStatus::IterationLimit && failure.iterations == 1,
            "Unconverged solve must report its iteration limit");
    const auto independent = matrix.multiply(failure.field);
    Vector error(rhs.size());
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        error[i] = rhs[i] - independent[i];
    }
    require(std::abs(failure.residual - vectorNorm(error) / vectorNorm(rhs)) < 1e-14,
            "Failure must report its true residual");
}

}  // namespace

int main() {
    try {
        testToleranceReport();
        const auto matrix = testMatrix();
#ifdef BPFEM_USE_MKL
        testMixedPrecisionGmres(matrix, expected);
#endif
        testIteration<GmresBackend>(matrix);
        testIteration<BiCGStabBackend>(matrix);
#ifdef BPFEM_USE_MKL
        const IterativeOperator sparseBlas(matrix, true);
        const auto rhs = matrix.multiply(expected);
        const auto actual = sparseBlas.multiply(expected);
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            require(std::abs(rhs[i] - actual[i]) < 1e-12, "MKL sparse operator changed matrix semantics");
        }
        GmresBackend solver;
        solver.setMklSpmv(true);
        solver.setPreconditioner(std::make_shared<PreconSchwarz>(
            2, 0, [] { return std::make_unique<PardisoBackend>(); }));
        check(solver.solve(matrix, rhs, strict), expected);
#endif
#ifdef BPFEM_USE_MUMPS
        for (const auto mode : {MumpsMode::InCore, MumpsMode::BlockLowRank, MumpsMode::OutOfCore}) {
            MumpsBackend mumps(mode, 1e-12);
            FactorizedSolveSession session(mumps, matrix, strict);
            const auto b = matrix.multiply(expected);
            const auto results = session.solveBatch({b, b});
            check(results.at(0), expected);
            check(results.at(1), expected);
            require(mumps.factorizationStatistics().numericFactorizationCount == 1,
                    "MUMPS batch did not reuse its factorization");
        }
        GmresBackend localMumps;
        localMumps.setPreconditioner(std::make_shared<PreconSchwarz>(
            2, 1, [] { return std::make_unique<MumpsBackend>(MumpsMode::InCore, 1e-12); }));
        check(localMumps.solve(matrix, matrix.multiply(expected), strict), expected);
#endif
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
