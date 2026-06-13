#pragma once

#include <filesystem>
#include <string>

namespace fem {

enum class SweepStrategy {
    Direct,  // legacy per-frequency direct solve
    Alps,    // Adaptive Lanczos-Pade Sweep (single-point Krylov MOR for now)
    Awe,     // Asymptotic Waveform Evaluation (single-point Padé MVP)
    Mgawe,   // Multipoint Galerkin AWE (global reduced basis)
    Wcawe    // Well-Conditioned AWE (orthogonalized moment basis)
};

enum class PortMethod {
    // Numerical single-mode 2D H(curl) eigensolve on the port face mesh.
    // Internally calls PortModeSolver::solve(faceId), which assembles the
    // port-face stiffness/mass and solves K_port v = k_c^2 M_port v. The
    // dominant (lowest k_c^2) physical eigenmode is used as the port mode.
    // Works for arbitrary cross sections and is the default.
    Numerical,
    // Analytic closed-form TE10 (rectangular only) projected onto FE port
    // DOFs (bpfem::apm). Kept as a legacy / reference path; aborts on
    // non-rectangular ports.
    Analytic,
    // Multi-mode transfinite element (bpfem::tfe). Same 2D eigensolve as
    // Numerical, but keeps the `tfeModesPerPort` lowest physical modes,
    // each contributing a rank-1 port operator.
    Transfinite
};

// Sparse-solver backend selection. Auto mirrors the legacy compile-time
// choice (PARDISO when built with BPFEM_USE_MKL, else BiCGSTAB) so unchanged
// CLI invocations keep their existing behavior. Direct / BiCGStab override
// the choice at runtime.
enum class LinearSolverKind {
    Auto,
    Direct,
    BiCGStab,
    Gmres
};

// Iterative-solver preconditioner selection. None = no preconditioning
// (default; matches behavior before Step 3 of strategy-interfaces phase).
// Jacobi = diagonal preconditioner.
// Ilu0 = zero-fill ILU(0).
// Direct backends ignore this field.
enum class PreconditionerKind {
    None,
    Jacobi,
    Ilu0
};

struct Options {
    std::filesystem::path aedt = "wg_bp_filter.aedt";
    std::filesystem::path mesh = "current.ngmesh";
    std::filesystem::path outDir = "results";
    int maxSweepPoints = 101;
    int maxIterations = 400;
    int basisOrder = 0;
    int fieldOutputOrder = 0;  // 0 means "auto: max(1, basisOrder + 1)"
    double tolerance = 1.0e-7;
    bool writeAllFields = true;

    // ALPS / fast-sweep options.
    SweepStrategy sweepStrategy = SweepStrategy::Direct;
    int alpsKrylovOrder = 30;             // q per port; ROM dim ~ Np * q after deflation
    double alpsExpansionFrequencyHz = 0.0;  // 0 means "use band center"
    int aweOrder = 8;                       // q for [q-1/q] Padé
    double aweExpansionFrequencyHz = 0.0;   // 0 means "use band center"
    int mgawePointCount = 3;                // number of expansion points
    int mgaweOrder = 8;                     // local moment count per expansion point
    double mgaweDropTolerance = 1.0e-10;    // global basis deflation threshold
    int wcaweOrder = 12;                    // target well-conditioned basis size
    double wcaweExpansionFrequencyHz = 0.0; // 0 means "use band center"
    double wcaweDropTolerance = 1.0e-12;    // R diagonal / MGS deflation threshold

    // Port-mode construction strategy. Default = Numerical (2D H(curl)
    // eigensolve on the port face mesh). Use Analytic only when an exact
    // TE10 closed-form reference is needed on a known rectangular port.
    PortMethod portMethod = PortMethod::Numerical;
    int tfeModesPerPort = 1;

    // Sparse solver backend selection (Step 1 of strategy-interfaces phase).
    LinearSolverKind linearSolver = LinearSolverKind::Auto;

    // Iterative-solver preconditioner selection (Step 3).
    PreconditionerKind preconditioner = PreconditionerKind::None;

    // GMRES restart length (only used when linearSolver == Gmres).
    int gmresRestart = 30;
};

Options parseOptions(int argc, char** argv);
int runApplication(int argc, char** argv);

}  // namespace fem
