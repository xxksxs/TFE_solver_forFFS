#include "bpfem/factory/SparseSolverFactory.hpp"
#include "bpfem/linalg/BiCGStabBackend.hpp"
#include "bpfem/linalg/FactorizedSolveSession.hpp"
#include "bpfem/linalg/GmresBackend.hpp"
#include "bpfem/linalg/IterativeOperator.hpp"
#include "bpfem/linalg/MumpsBackend.hpp"
#include "bpfem/linalg/PardisoBackend.hpp"
#include "bpfem/linalg/PreconSchwarz.hpp"

#include <cmath>
#include <future>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>

using namespace fem;
using namespace fem::linalg;

namespace {

using Vector = std::vector<std::complex<double>>;

static_assert(!std::is_copy_constructible_v<GmresBackend> && std::is_move_constructible_v<GmresBackend>);
static_assert(!std::is_copy_constructible_v<BiCGStabBackend> && std::is_move_constructible_v<BiCGStabBackend>);

class CountingIdentity final : public IPreconditioner {
public:
    void setup(const SparseMatrix&) override {
        if (failSetup) {
            throw std::runtime_error("Injected setup failure");
        }
        ++preparations;
    }
    void apply(const Vector& x, Vector& y) const override {
        y = x;
        if (nonFiniteOutput) {
            y.front() = std::numeric_limits<double>::quiet_NaN();
        }
    }
    void clear() noexcept override { ++clears; }
    const char* name() const override { return "counting identity"; }
    int preparations = 0;
    int clears = 0;
    bool failSetup = false;
    bool nonFiniteOutput = false;
};

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

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
const Vector otherExpected{{-0.8, 0.0}, {0.2, 0.8}, {0.0, -1.0}, {0.4, -0.5}};
const SolverConfig strict{80, 1e-11};

void check(const SolveResult& result, const Vector& solution) {
    require(result.status == SolveStatus::Converged, "Solve did not converge");
    require(result.residual <= strict.tolerance, "Original-matrix residual is too large");
    require(result.field.size() == solution.size(), "Incorrect solution size");
    for (std::size_t i = 0; i < solution.size(); ++i) {
        require(std::abs(result.field[i] - solution[i]) < 1e-9, "Wrong complex solution");
    }
}

void checkSessionIsolation(ISparseSolver& solver, const SparseMatrix& matrix) {
    auto next = matrix;
    next.add(0, 0, {0.5, 0.9});
    {
        FactorizedSolveSession first(solver, matrix, strict);
        check(first.solve(matrix.multiply(expected)), expected);
        bool rejected = false;
        try {
            FactorizedSolveSession overlapping(solver, next, strict);
        } catch (const std::logic_error&) {
            rejected = true;
        }
        require(rejected, "A second session can replace the active matrix");
        check(first.solve(matrix.multiply(otherExpected)), otherExpected);
    }
    FactorizedSolveSession second(solver, next, strict);
    check(second.solve(next.multiply(expected)), expected);
}

void checkConcurrentSessionIsolation(ISparseSolver& solver, const SparseMatrix& matrix) {
    auto next = matrix;
    next.add(0, 0, {0.5, 0.9});
    const auto rhs = matrix.multiply(expected);
    const auto nextRhs = next.multiply(otherExpected);
    std::promise<void> start;
    std::promise<void> localAttempted;
    auto startSignal = start.get_future();
    auto localAttemptedSignal = localAttempted.get_future();
    auto open = [&](const SparseMatrix& coefficients) {
        try {
            return std::make_unique<FactorizedSolveSession>(solver, coefficients, strict);
        } catch (const std::logic_error&) {
            return std::unique_ptr<FactorizedSolveSession>{};
        }
    };
    auto worker = std::async(std::launch::async, [&] {
        startSignal.wait();
        auto session = open(next);
        localAttemptedSignal.wait();
        if (session) {
            check(session->solve(nextRhs), otherExpected);
        }
        return static_cast<bool>(session);
    });

    start.set_value();
    std::unique_ptr<FactorizedSolveSession> local;
    try {
        local = open(matrix);
    } catch (...) {
        localAttempted.set_value();
        throw;
    }
    localAttempted.set_value();
    const bool workerWon = worker.get();
    require(static_cast<bool>(local) != workerWon, "Exactly one concurrent session must acquire the solver");
    if (local) {
        check(local->solve(rhs), expected);
        local.reset();
    }
    FactorizedSolveSession reopened(solver, next, strict);
    check(reopened.solve(nextRhs), otherExpected);
}

void checkFactoryPreconditionerIsolation(LinearSolverKind kind) {
    Options options;
    options.linearSolver = kind;
    options.preconditioner = PreconditionerKind::Jacobi;
    auto firstSolver = factory::makeSparseSolver(options);
    auto secondSolver = factory::makeSparseSolver(options);
    SparseMatrix firstMatrix(4), secondMatrix(4);
    for (int i = 0; i < 4; ++i) {
        firstMatrix.add(i, i, {2.0 + static_cast<double>(i), 0.4});
        secondMatrix.add(i, i, {3.0, 1.0 + static_cast<double>(i)});
    }
    const SolverConfig oneIteration{1, strict.tolerance};
    FactorizedSolveSession first(*firstSolver, firstMatrix, oneIteration);
    {
        FactorizedSolveSession second(*secondSolver, secondMatrix, oneIteration);
        check(first.solve(firstMatrix.multiply(expected)), expected);
        check(second.solve(secondMatrix.multiply(otherExpected)), otherExpected);
    }
    check(first.solve(firstMatrix.multiply(otherExpected)), otherExpected);
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
        const auto results = session.solveBatch({rhs, zero, matrix.multiply(otherExpected)});
        check(results.at(0), expected);
        check(results.at(1), zero);
        check(results.at(2), otherExpected);
        constexpr double rhsScale = 1e-24;
        auto smallRhs = rhs;
        for (auto& entry : smallRhs) {
            entry *= rhsScale;
        }
        auto smallResult = session.solve(smallRhs);
        for (auto& entry : smallResult.field) {
            entry /= rhsScale;
        }
        check(smallResult, expected);
        require(results[1].iterations == 0, "Zero RHS should require no iterations");
        require(precon->preparations == 1, "Prepared solve repeated preconditioner setup");
    }
    require(precon->clears == 1, "Session exit did not clear the prepared preconditioner");
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
    const SparseMatrix singular(matrix.size());
    const auto breakdown = solver.solve(singular, rhs, strict);
    require(breakdown.status == SolveStatus::Breakdown && breakdown.residual == 1.0,
            "Inconsistent zero matrix must report breakdown with unit true residual");
    auto invalid = matrix;
    invalid.add(0, 0, {std::numeric_limits<double>::quiet_NaN(), 0.0});
    const auto nonFinite = solver.solve(invalid, rhs, strict);
    require(nonFinite.status == SolveStatus::NonFinite && !std::isfinite(nonFinite.residual),
            "Non-finite matrix must report a non-finite true residual");
    require(solver.solve(invalid, zero, strict).status == SolveStatus::NonFinite,
            "Zero RHS must also be checked on the original matrix");
    precon->nonFiniteOutput = true;
    const auto invalidPreconditioner = solver.solve(matrix, rhs, strict);
    require(invalidPreconditioner.status == SolveStatus::NonFinite
                && invalidPreconditioner.residual == 1.0,
            "Non-finite preconditioner must be reported even when the last iterate is finite");
    precon->nonFiniteOutput = false;
    precon->failSetup = true;
    bool rejected = false;
    try {
        FactorizedSolveSession invalidSession(solver, matrix, strict);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    require(rejected, "Session must propagate preparation errors");
    precon->failSetup = false;
    checkSessionIsolation(solver, matrix);
}

}  // namespace

int main() {
    try {
        const auto matrix = testMatrix();
        testIteration<GmresBackend>(matrix);
        testIteration<BiCGStabBackend>(matrix);
        GmresBackend concurrent;
        checkConcurrentSessionIsolation(concurrent, matrix);
        checkFactoryPreconditionerIsolation(LinearSolverKind::Gmres);
        checkFactoryPreconditionerIsolation(LinearSolverKind::BiCGStab);
#ifdef BPFEM_USE_MKL
        PardisoBackend pardiso;
        checkSessionIsolation(pardiso, matrix);
        checkConcurrentSessionIsolation(pardiso, matrix);
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
            checkSessionIsolation(mumps, matrix);
            const auto factorizationsBefore = mumps.factorizationStatistics().numericFactorizationCount;
            FactorizedSolveSession session(mumps, matrix, strict);
            const auto b = matrix.multiply(expected);
            const auto results = session.solveBatch({b, matrix.multiply(otherExpected)});
            check(results.at(0), expected);
            check(results.at(1), otherExpected);
            require(mumps.factorizationStatistics().numericFactorizationCount == factorizationsBefore + 1,
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
