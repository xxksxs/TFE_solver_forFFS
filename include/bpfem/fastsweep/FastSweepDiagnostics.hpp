#pragma once

#include "bpfem/core/Types.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace fem::fastsweep {

struct FastSweepDiagnostics {
    std::string algorithm;
    std::vector<double> expansionFrequenciesHz;
    int requestedOrder = 0;
    int romDimension = 0;
    int retainedColumns = 0;
    int deflatedColumns = 0;
    double basisOrthogonalityError = 0.0;
    double padeInputPivotRatio = 0.0;
    double padeOutputPivotRatio = 0.0;
    double wcaweMomentReconstructionError = 0.0;
    bool reducedSolveSucceeded = true;
    double maxPassivityError = 0.0;
};

double maxPassivityError(const std::vector<SParameterPoint>& points);

bool writeDiagnosticsJson(const std::filesystem::path& path,
                          const FastSweepDiagnostics& diagnostics);

}  // namespace fem::fastsweep
