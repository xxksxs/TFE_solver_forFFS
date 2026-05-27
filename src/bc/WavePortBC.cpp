#include "bpfem/bc/WavePortBC.hpp"

#include "bpfem/core/Constants.hpp"

#include <algorithm>
#include <cmath>
#include <complex>

namespace fem::bc {

namespace {

// Apply one analytic mode of a port: rank-1 update on `matrix` plus optional
// incident-amplitude contribution on `rhs`. Mirrors the legacy
// FEMAssembler::applyWavePorts inner lambda 1:1 so output stays bit-for-bit
// identical for the existing default path.
void applyOneMode(SparseMatrix& matrix,
                  std::vector<std::complex<double>>& rhs,
                  const PortMode& portMode,
                  bool excited,
                  double magnitudeW,
                  double phaseDeg,
                  double frequencyHz,
                  const std::unordered_set<int>& constrainedDofs) {
    const auto& mode = portMode.couplingWeights;
    if (mode.empty()) {
        return;
    }
    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double s = powerNormalizationFactor(portMode, frequencyHz);
    const std::complex<double> incident = (excited && s > 0.0)
        ? std::polar(std::sqrt(std::max(magnitudeW, 0.0)) * s, phaseDeg * pi / 180.0)
        : std::complex<double>(0.0, 0.0);
    const double beta = std::sqrt(std::max(0.0, k0 * k0 - portMode.cutoffWavenumberSquared));
    const std::complex<double> admittance(0.0, beta > 0.0 ? beta : k0);
    for (std::size_t rowIndex = 0; rowIndex < mode.size(); ++rowIndex) {
        const auto& [row, rowValue] = mode[rowIndex];
        if (constrainedDofs.count(row) != 0) {
            continue;
        }
        rhs[static_cast<std::size_t>(row)] += 2.0 * admittance * incident * rowValue;
        for (std::size_t colIndex = rowIndex; colIndex < mode.size(); ++colIndex) {
            const auto& [col, colValue] = mode[colIndex];
            if (constrainedDofs.count(col) != 0) {
                continue;
            }
            matrix.add(row, col, admittance * rowValue * colValue);
        }
    }
}

}  // namespace

void WavePortBC::apply(SparseMatrix& matrix,
                       std::vector<std::complex<double>>& rhs,
                       const AssemblyContext& ctx) const {
    for (const auto& port : ctx.project.ports) {
        // TFE multi-mode path: each analytic mode contributes a rank-1 port
        // operator. Only the excitation mode of an excited project port
        // carries an incident amplitude; the other modes are pure absorbers.
        if (const auto* multi = ctx.portModeSolver.multiMode(port.faceId); multi != nullptr) {
            for (std::size_t mIdx = 0; mIdx < multi->modes.size(); ++mIdx) {
                const bool isExcited = port.excited
                    && (static_cast<int>(mIdx) == multi->excitationModeIndex);
                applyOneMode(matrix, rhs, multi->modes[mIdx], isExcited,
                             port.magnitudeW, port.phaseDeg, ctx.frequencyHz,
                             ctx.constrainedDofs);
            }
            continue;
        }
        // Single-mode path (NPM, APM): unchanged behavior.
        applyOneMode(matrix, rhs, ctx.portModeSolver.solve(port.faceId),
                     port.excited, port.magnitudeW, port.phaseDeg,
                     ctx.frequencyHz, ctx.constrainedDofs);
    }
}

void WavePortBC::declareSparsity(SparsePatternBuilder& builder,
                                 const AssemblyContext& ctx) const {
    auto declareCoupling = [&](const std::vector<std::pair<int, double>>& mode) {
        for (std::size_t rowIndex = 0; rowIndex < mode.size(); ++rowIndex) {
            const int row = mode[rowIndex].first;
            for (std::size_t colIndex = rowIndex; colIndex < mode.size(); ++colIndex) {
                builder.add(row, mode[colIndex].first);
            }
        }
    };
    for (const auto& port : ctx.project.ports) {
        if (const auto* multi = ctx.portModeSolver.multiMode(port.faceId); multi != nullptr) {
            for (const auto& m : multi->modes) {
                declareCoupling(m.couplingWeights);
            }
        } else {
            declareCoupling(ctx.portModeSolver.solve(port.faceId).couplingWeights);
        }
    }
}

std::vector<AffinePortContribution> WavePortBC::affineContributions(
    const AssemblyContext& ctx) const {
    std::vector<AffinePortContribution> out;
    out.reserve(ctx.project.ports.size());
    for (std::size_t projectIdx = 0; projectIdx < ctx.project.ports.size(); ++projectIdx) {
        const auto& port = ctx.project.ports[projectIdx];
        auto pushOneMode = [&](const PortMode& portMode, bool excited) {
            AffinePortContribution c;
            c.coupling.reserve(portMode.couplingWeights.size());
            for (const auto& [dof, weight] : portMode.couplingWeights) {
                if (ctx.constrainedDofs.count(dof) != 0) {
                    continue;
                }
                c.coupling.emplace_back(dof, weight);
            }
            c.cutoffSquared = portMode.cutoffWavenumberSquared;
            c.faceId = port.faceId;
            c.projectPortIndex = static_cast<int>(projectIdx);
            c.isExcitationMode = excited;
            out.push_back(std::move(c));
        };

        if (const auto* multi = ctx.portModeSolver.multiMode(port.faceId); multi != nullptr) {
            for (std::size_t mIdx = 0; mIdx < multi->modes.size(); ++mIdx) {
                const bool isExc = port.excited
                    && (static_cast<int>(mIdx) == multi->excitationModeIndex);
                pushOneMode(multi->modes[mIdx], isExc);
            }
        } else {
            pushOneMode(ctx.portModeSolver.solve(port.faceId), port.excited);
        }
    }
    return out;
}

}  // namespace fem::bc

