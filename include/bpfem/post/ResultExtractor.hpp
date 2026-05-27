#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fem/PortModeSolver.hpp"

#include <complex>
#include <vector>

namespace fem {

class ResultExtractor {
public:
    ResultExtractor(const ProjectDefinition& project, const PortModeSolver& portModeSolver);

    SParameterPoint extract(double frequencyHz, const std::vector<std::complex<double>>& edgeDofs) const;

private:
    const ProjectDefinition& project_;
    const PortModeSolver& portModeSolver_;

    std::complex<double> portProjection(const std::vector<std::complex<double>>& edgeDofs, int faceId) const;
};

}  // namespace fem
