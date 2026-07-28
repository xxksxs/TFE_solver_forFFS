#pragma once

#include "bpfem/core/Logger.hpp"
#include "bpfem/core/Types.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/post/ResultExtractor.hpp"

#include <complex>
#include <filesystem>
#include <functional>
#include <vector>

namespace fem::sweep {

// Per-call data passed to a sweep strategy. Holds references only.
//
// solver 是应用层选定的线性求解后端。ALPS 与 AWE-family 通过能力接口
// 在离线阶段复用展开点分解，迭代后端则自动回退到普通 solve。
struct SweepContext {
    const ProjectDefinition& project;
    const FEMAssembler& assembler;
    const PortModeSolver& portModeSolver;
    const ResultExtractor& extractor;
    linalg::ISparseSolver& solver;
    Logger& log;
    std::filesystem::path outputDirectory;

    // Per-frequency callback invoked by direct-style strategies after a
    // successful solve. Fast-sweep ROM strategies usually ignore it and only
    // publish SweepResult::lastEdgeDofs for field_last.vtu, avoiding one
    // reconstructed full-space field per frequency. Default empty closure =
    // no field output.
    using FieldSink = std::function<void(double frequencyHz,
                                         const std::vector<std::complex<double>>& edgeDofs)>;
    FieldSink onFieldSolved;

    // Per-frequency progress callback. DirectSweep invokes this after each
    // assemble + solve with timing breakdown so Application can publish
    // [stat] log rows and append to the JSON sweep array. AlpsSweep does
    // not call it (its online evaluate is too cheap to merit per-point
    // breakdown). Default empty closure = no progress instrumentation.
    using ProgressSink = std::function<void(double frequencyHz,
                                            double assembleSec,
                                            double solveSec,
                                            int iterations,
                                            double residual)>;
    ProgressSink onFrequencyStat;

    // Linear-solver tolerance / iteration budget forwarded by Application.
    // Direct backends ignore both.
    int linearMaxIterations = 400;
    double linearTolerance = 1.0e-7;
};

// Output of a sweep run.
//
// `lastEdgeDofs` is populated when a strategy can provide a representative
// final full-space field so Application can write `field_last.vtu`.
struct SweepResult {
    std::vector<SParameterPoint> points;
    std::vector<std::complex<double>> lastEdgeDofs;
    double lastFrequencyHz = 0.0;

    // Fast-sweep 的细粒度性能统计；direct 等策略可保持为零。
    double offlineBuildSec = 0.0;
    double portLinearizationSec = 0.0;
    double lanczosOperatorSec = 0.0;
    double poleDecompositionSec = 0.0;
    double orthogonalizationSec = 0.0;
    double romProjectionSec = 0.0;
    double onlineSweepSec = 0.0;
};

// Frequency-sweep strategy interface.
//
// Concrete kinds:
//   - DirectSweep   per-frequency assemble + ISparseSolver::solve, no MOR
//   - AlpsSweep     single-point two-sided Lanczos-Pade transfer models
class ISweepStrategy {
public:
    virtual ~ISweepStrategy() = default;

    // `frequencies` is the user-supplied list. Implementations may evaluate
    // those points exactly (DirectSweep) or use them as a basis for adaptive
    // refinement and return a different set (future AdaptiveSweep). The
    // returned `SweepResult::points` must be sorted by frequency.
    virtual SweepResult run(const std::vector<double>& frequencies,
                            const SweepContext& ctx) = 0;

    virtual const char* name() const = 0;
};

}  // namespace fem::sweep

