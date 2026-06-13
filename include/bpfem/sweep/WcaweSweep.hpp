#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"
#include "bpfem/sweep/ISweepStrategy.hpp"

#include <complex>
#include <filesystem>
#include <vector>

namespace fem::sweep {

struct WcaweOptions {
    // Target basis size. The builder may stop earlier if MGS detects a
    // near-dependent moment direction.
    int order = 12;

    // Expansion frequency in Hz. 0 means "use band center".
    double expansionFrequencyHz = 0.0;

    // Deflation threshold for the MGS residual / R diagonal.
    double dropTolerance = 1.0e-12;
};

class WcaweSweep : public ISweepStrategy {
public:
    using Complex = std::complex<double>;

    WcaweSweep(const ProjectDefinition& project,
               const FEMAssembler& assembler,
               const PortModeSolver& portModeSolver,
               WcaweOptions options = {});

    SweepResult run(const std::vector<double>& frequencies,
                    const SweepContext& ctx) override;

    const char* name() const override { return "wcawe"; }

    int buildOffline(double expansionFrequencyHz,
                     linalg::ISparseSolver& solver,
                     const linalg::SolverConfig& solverConfig);
    SParameterPoint evaluate(double frequencyHz) const;
    std::vector<Complex> reconstructField(double frequencyHz) const;

    int dimension() const { return romDim_; }
    int deflatedColumns() const { return deflatedColumns_; }
    bool ready() const { return ready_; }

    bool writeBasisConditionCsv(const std::filesystem::path& path) const;

private:
    struct ConditionRecord {
        int order = 0;
        double aweConditionProxy = 0.0;
        double wcaweConditionProxy = 1.0;
        double rDiagonalAbs = 0.0;
        double orthogonalityError = 0.0;
    };

    std::vector<std::vector<Complex>> buildAweMoments(
        double expansionFrequencyHz,
        linalg::ISparseSolver& solver,
        const linalg::SolverConfig& solverConfig,
        const std::vector<std::vector<Complex>>& portVectors) const;

    bool appendWellConditionedColumn(const std::vector<Complex>& moment);
    void projectReducedModel(const std::vector<std::vector<Complex>>& portVectors);
    std::vector<Complex> solveReduced(double frequencyHz) const;
    const PortMode& virtualPortMode(int virtualPortIndex) const;

    const ProjectDefinition& project_;
    const FEMAssembler& assembler_;
    const PortModeSolver& portModeSolver_;
    WcaweOptions options_;

    bool ready_ = false;
    int romDim_ = 0;
    int fullDim_ = 0;
    int numProjectPorts_ = 0;
    int numVirtualPorts_ = 0;
    int deflatedColumns_ = 0;
    double expansionFrequencyHz_ = 0.0;

    FEMAssembler::AffineSystem affine_;
    std::vector<std::vector<Complex>> basis_;
    std::vector<Complex> triangularR_;
    std::vector<double> rDiagonalAbs_;
    std::vector<ConditionRecord> conditionRecords_;
    std::vector<Complex> Ktilde_;
    std::vector<Complex> Mtilde_;
    std::vector<std::vector<Complex>> portModeReduced_;
};

}  // namespace fem::sweep
