#include "bpfem/post/ResultExtractor.hpp"

#include "bpfem/core/Constants.hpp"

#include <cmath>

namespace fem {

ResultExtractor::ResultExtractor(const ProjectDefinition& project, const PortModeSolver& portModeSolver)
    : project_(project), portModeSolver_(portModeSolver) {}

SParameterPoint ResultExtractor::extract(double frequencyHz, const std::vector<std::complex<double>>& edgeDofs) const {
    if (project_.ports.size() < 2) {
        SParameterPoint sp;
        sp.frequencyHz = frequencyHz;
        return sp;
    }
    return extractFromPortProjections(frequencyHz,
        portProjection(edgeDofs, project_.ports[0].faceId),
        portProjection(edgeDofs, project_.ports[1].faceId));
}

SParameterPoint ResultExtractor::extractFromPortProjections(double frequencyHz,
                                                           std::complex<double> inputProjection,
                                                           std::complex<double> outputProjection) const {
    SParameterPoint sp;
    sp.frequencyHz = frequencyHz;
    if (project_.ports.size() < 2) {
        sp.s11 = 0.0;
        sp.s21 = 0.0;
        return sp;
    }
    const auto& mIn = portModeSolver_.solve(project_.ports[0].faceId);
    const auto& mOut = portModeSolver_.solve(project_.ports[1].faceId);
    const double sIn = powerNormalizationFactor(mIn, frequencyHz);
    const double sOut = powerNormalizationFactor(mOut, frequencyHz);
    if (sIn <= 0.0) {
        sp.s11 = 0.0;
        sp.s21 = 0.0;
        return sp;
    }
    const std::complex<double> bIn = inputProjection / sIn;
    const std::complex<double> bOut = (sOut > 0.0)
        ? outputProjection / sOut
        : std::complex<double>(0.0, 0.0);
    const std::complex<double> incident = project_.ports[0].excited
        ? std::polar(std::sqrt(std::max(project_.ports[0].magnitudeW, 0.0)), project_.ports[0].phaseDeg * pi / 180.0)
        : std::complex<double>(1.0, 0.0);
    sp.s11 = incident == std::complex<double>(0.0, 0.0) ? 0.0 : (bIn - incident) / incident;
    sp.s21 = incident == std::complex<double>(0.0, 0.0) ? 0.0 : bOut / incident;
    return sp;
}

std::complex<double> ResultExtractor::portProjection(const std::vector<std::complex<double>>& edgeDofs, int faceId) const {
    const auto& mode = portModeSolver_.solve(faceId).couplingWeights;
    if (mode.empty()) {
        return 0.0;
    }
    std::complex<double> sum = 0.0;
    for (const auto& [edgeIndex, value] : mode) {
        if (edgeIndex >= 0 && static_cast<std::size_t>(edgeIndex) < edgeDofs.size()) {
            sum += edgeDofs[static_cast<std::size_t>(edgeIndex)] * value;
        }
    }
    return sum;
}

}  // namespace fem
