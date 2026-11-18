#include "bpfem/sweep/MgaweSweep.hpp"

#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"
#include "bpfem/fastsweep/PortModeUtilities.hpp"

#include <algorithm>
#include <chrono>
#include <complex>
#include <cmath>
#include <limits>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fem::sweep {

namespace {
constexpr double localSpaceTolerance = 1.0e-6;
constexpr double residualTolerance = 1.0e-2;
using Vector = std::vector<std::complex<double>>;

// Slone et al. (2002), equations (20)-(24).
bool localSpaceExhausted(const SparseMatrix& a0, const std::vector<Vector>& basis,
                         const Vector& rhs) {
    const double rhsNorm = fastsweep::norm2(rhs);
    if (rhsNorm == 0.0) return true;
    const int q = static_cast<int>(basis.size());
    std::vector<std::complex<double>> matrix(static_cast<std::size_t>(q) * q);
    Vector coordinates(static_cast<std::size_t>(q));
    std::vector<Vector> applied;
    applied.reserve(basis.size());
    for (const auto& v : basis) applied.push_back(a0.multiply(v));
    for (int i = 0; i < q; ++i) {
        coordinates[i] = fastsweep::bilinear(basis[i], rhs);
        for (int j = 0; j < q; ++j) {
            matrix[static_cast<std::size_t>(i) * q + j] = fastsweep::bilinear(basis[i], applied[j]);
        }
    }
    if (!fastsweep::denseSolve(matrix, coordinates, q)) return false;
    auto remainder = rhs;
    for (int j = 0; j < q; ++j) {
        for (std::size_t i = 0; i < rhs.size(); ++i) remainder[i] -= applied[j][i] * coordinates[j];
    }
    return fastsweep::norm2(remainder) <= localSpaceTolerance * rhsNorm;
}

struct ExpansionChoice {
    double frequency = 0.0;
    double maxResidual = 0.0;
};

// Slone et al. (2002), equations (19) and (25)-(27), on the requested grid.
ExpansionChoice nextExpansion(const fastsweep::GalerkinReducedModel& model,
                              const std::vector<double>& frequencies,
                              const std::vector<double>& expansions) {
    const int count = static_cast<int>(frequencies.size());
    std::vector<double> residuals;
    residuals.reserve(frequencies.size());
    for (double f : frequencies) residuals.push_back(model.relativeResidual(f));
    const auto worst = std::max_element(residuals.begin(), residuals.end());
    const double maxResidual = *worst;
    if (maxResidual < residualTolerance) return {0.0, maxResidual};

    auto used = [&](int i) {
        return std::find(expansions.begin(), expansions.end(), frequencies[i]) != expansions.end();
    };
    int peak = -1;
    for (int i = 0; i < count; ++i) {
        if (!used(i) && residuals[i] >= residualTolerance
            && (peak < 0 || residuals[i] > residuals[peak])) peak = i;
    }
    if (peak < 0) return {0.0, maxResidual};
    int left = peak;
    int right = peak;
    while (left > 0 && residuals[left - 1] >= residualTolerance && !used(left - 1)) --left;
    while (right + 1 < count && residuals[right + 1] >= residualTolerance && !used(right + 1)) ++right;

    const bool touchesLeft = left == 0;
    const bool touchesRight = right == count - 1;
    const int leftBoundary = touchesLeft ? -1 : left - 1;
    const int rightBoundary = touchesRight ? count : right + 1;
    const double leftWeight = touchesLeft == touchesRight ? 0.5 : (touchesLeft ? 0.75 : 0.25);
    const int next = std::clamp(static_cast<int>(std::lround(
        leftWeight * leftBoundary + (1.0 - leftWeight) * rightBoundary)), left, right);
    return {frequencies[next], maxResidual};
}
}  // namespace

MgaweSweep::MgaweSweep(const ProjectDefinition& project,
                       const FEMAssembler& assembler,
                       const PortModeSolver& portModeSolver,
                       MgaweOptions options)
    : project_(project),
      assembler_(assembler),
      portModeSolver_(portModeSolver),
      options_(options) {
    if (options_.adaptive && !options_.expansions.empty()) {
        throw std::invalid_argument("MGAWE: explicit expansions are only supported in fixed-point modes");
    }
    if (options_.adaptive && options_.useLastPolynomial) {
        throw std::invalid_argument("MGAWE: --mgawe-2001 and --mgawe-adaptive are mutually exclusive");
    }
    options_.expansionPointCount = std::max(1, options_.expansionPointCount);
    options_.localOrder = std::max(1, options_.localOrder);
    if (options_.dropTolerance <= 0.0) {
        options_.dropTolerance = 1.0e-10;
    }
}

SweepResult MgaweSweep::run(const std::vector<double>& frequencies, const SweepContext& ctx) {
    if (frequencies.empty()) {
        return {};
    }

    const double left = frequencies.front();
    const double right = frequencies.back();
    const int pointCount = std::max(1, options_.expansionPointCount);
    std::vector<double> expansionFrequencies;
    expansionFrequencies.reserve(static_cast<std::size_t>(pointCount));
    if (pointCount == 1 || left == right) {
        expansionFrequencies.push_back(0.5 * (left + right));
    } else {
        for (int i = 0; i < pointCount; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(pointCount - 1);
            expansionFrequencies.push_back(left + t * (right - left));
        }
    }

    if (!options_.expansions.empty()) {
        expansionFrequencies.clear();
        for (const auto& [frequency, order] : options_.expansions) {
            if (!std::isfinite(frequency) || frequency <= 0.0 || order < 1
                || std::find(expansionFrequencies.begin(), expansionFrequencies.end(), frequency)
                    != expansionFrequencies.end()) {
                throw std::invalid_argument("MGAWE explicit expansion points must be positive and distinct");
            }
            expansionFrequencies.push_back(frequency);
        }
    }

    std::ostringstream points;
    points << std::fixed << std::setprecision(6);
    for (std::size_t i = 0; i < expansionFrequencies.size(); ++i) {
        if (i != 0) {
            points << ", ";
        }
        points << (expansionFrequencies[i] / 1.0e9);
    }
    if (!options_.adaptive) ctx.log.info("MGAWE expansion frequencies: [" + points.str() + "] GHz");
    ctx.log.info("MGAWE local AWE moment count per point: " + std::to_string(options_.localOrder));
    std::ostringstream dropTol;
    dropTol << std::scientific << std::setprecision(3) << options_.dropTolerance;
    ctx.log.info("MGAWE global MGS drop tolerance: " + dropTol.str());

    linalg::SolverConfig solverCfg;
    solverCfg.maxIterations = ctx.linearMaxIterations;
    solverCfg.tolerance = ctx.linearTolerance;
    const auto offlineStarted = std::chrono::steady_clock::now();
    if (options_.adaptive) {
        expansionFrequencies = {frequencies[(frequencies.size() - 1) / 2]};
    }
    const int dim = build(expansionFrequencies, frequencies, ctx.solver, solverCfg);
    if (options_.adaptive) {
        std::ostringstream summary;
        summary << "MGAWE adaptive stop: "
                << (adaptiveConverged_ ? "residual tolerance" : "configured expansion-point limit")
                << ", max relative residual " << maxRelativeResidual_ << ", expansions:";
        for (std::size_t i = 0; i < expansionFrequencies.size(); ++i) {
            summary << ' ' << expansionFrequencies[i] / 1.0e9 << "GHz/q=" << localOrders_[i];
        }
        ctx.log.info(summary.str());
    }
    const double offlineBuildSec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - offlineStarted).count();
    ctx.log.info("MGAWE ROM dimension: " + std::to_string(dim));
    ctx.log.info("MGAWE retained candidate columns: " + std::to_string(retainedColumns_));
    ctx.log.info("MGAWE deflated candidate columns: " + std::to_string(deflatedColumns_));

    SweepResult out;
    out.offlineBuildSec = offlineBuildSec;
    out.orthogonalizationSec = orthogonalizationSec_;
    out.romProjectionSec = romProjectionSec_;
    out.points.reserve(frequencies.size());
    const auto onlineStarted = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        const double f = frequencies[i];
        if (i == 0 || (i + 1) % 50 == 0 || i + 1 == frequencies.size()) {
            ctx.log.info("MGAWE evaluating point " + std::to_string(i + 1) + "/"
                         + std::to_string(frequencies.size())
                         + " @ " + std::to_string(f / 1.0e9) + " GHz");
        }
        out.points.push_back(evaluate(f));
        out.lastFrequencyHz = f;
    }
    out.onlineSweepSec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - onlineStarted).count();
    if (!frequencies.empty()) {
        out.lastEdgeDofs = reconstructField(frequencies.back());
    }
    if (!ctx.outputDirectory.empty()) {
        fastsweep::FastSweepDiagnostics diag;
        diag.algorithm = name();
        diag.model = options_.polynomialProjection ? "wavenumber_taylor_polynomial" : "exact_port_frequency_law";
        if (options_.useLastPolynomial) diag.model = "wavenumber_taylor_last_expansion_2001";
        diag.localOrders = localOrders_;
        diag.stoppingReason = options_.adaptive
            ? (adaptiveConverged_ ? "residual_tolerance" : "expansion_budget_or_no_unused_point")
            : "fixed_expansion_points";
        if (options_.adaptive) diag.maximumRelativeResidual = maxRelativeResidual_;
        diag.expansionFrequenciesHz = expansionFrequencies;
        diag.requestedOrder = options_.localOrder;
        diag.romDimension = romDim_;
        diag.retainedColumns = retainedColumns_;
        diag.deflatedColumns = deflatedColumns_;
        diag.basisOrthogonalityError = model_.basisOrthogonalityError();
        diag.reducedSolveSucceeded = model_.lastSolveSucceeded();
        diag.maxPassivityError = fastsweep::maxPassivityError(out.points);
        const auto path = ctx.outputDirectory / "diagnostics.json";
        if (fastsweep::writeDiagnosticsJson(path, diag)) {
            ctx.log.info("Wrote " + path.string());
        } else {
            ctx.log.warn("Failed to write " + path.string());
        }
    }
    (void)ctx.onFieldSolved;
    return out;
}

int MgaweSweep::buildOffline(const std::vector<double>& expansionFrequencies,
                             linalg::ISparseSolver& solver,
                             const linalg::SolverConfig& solverConfig) {
    auto points = expansionFrequencies;
    return build(points, {}, solver, solverConfig);
}

int MgaweSweep::build(std::vector<double>& expansionFrequencies,
                      const std::vector<double>& sweepFrequencies,
                      linalg::ISparseSolver& solver, const linalg::SolverConfig& solverConfig) {
    ready_ = false;
    polynomials_.clear();
    localOrders_.clear();
    adaptiveConverged_ = false;
    maxRelativeResidual_ = 0.0;
    if (project_.ports.size() < 2) {
        throw std::runtime_error("MgaweSweep: at least two project ports are required");
    }
    if (expansionFrequencies.empty()) {
        throw std::runtime_error("MgaweSweep: at least one expansion point is required");
    }

    numProjectPorts_ = static_cast<int>(project_.ports.size());
    affine_ = assembler_.buildAffineSystem();
    fullDim_ = static_cast<int>(affine_.K.size());
    numVirtualPorts_ = static_cast<int>(affine_.portCoupling.size());
    if (numVirtualPorts_ == 0) {
        throw std::runtime_error("MgaweSweep: AffineSystem produced zero port entries");
    }
    if (!affine_.lossless) {
        throw std::runtime_error("MgaweSweep: lossy materials (sigma != 0) are not supported. Use --sweep direct.");
    }

    const auto portVectors =
        fastsweep::PolynomialPortMomentBuilder::buildPortVectors(affine_, affine_.K.size());

    basis_.clear();
    retainedColumns_ = 0;
    deflatedColumns_ = 0;
    orthogonalizationSec_ = 0.0;
    romProjectionSec_ = 0.0;
    const int candidateCount =
        static_cast<int>(expansionFrequencies.size()) * std::max(1, options_.localOrder);
    basis_.reserve(static_cast<std::size_t>(candidateCount));

    for (std::size_t point = 0; point < expansionFrequencies.size(); ++point) {
        const double expansionHz = expansionFrequencies[point];
        const bool adaptive = options_.adaptive && !sweepFrequencies.empty();
        using Builder = fastsweep::PolynomialPortMomentBuilder;
        const auto variable = options_.polynomialProjection ? Builder::Variable::WavenumberShift
                                                            : Builder::Variable::LambdaShift;
        const int requestedOrder = !options_.expansions.empty()
            ? options_.expansions.at(point).second : options_.localOrder;
        const auto moments = Builder::generateLosslessMoments(
            project_, assembler_, portModeSolver_, affine_, expansionHz, requestedOrder,
            solver, solverConfig, portVectors, variable,
            [&](const SparseMatrix& a0, const Vector& rhs) {
                return !adaptive || basis_.empty() || !localSpaceExhausted(a0, basis_, rhs);
            },
            [&](const Vector& moment) {
                const auto start = std::chrono::steady_clock::now();
                appendIfIndependent(Vector(moment));
                orthogonalizationSec_ += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start).count();
            });
        localOrders_.push_back(static_cast<int>(moments.size()));
        if (options_.polynomialProjection) {
            polynomials_.push_back(Builder::buildTaylorSeries(project_, portModeSolver_, affine_,
                expansionHz, std::max(3, requestedOrder), Builder::Variable::WavenumberShift));
        }
        if (adaptive && !basis_.empty()) {
            const auto start = std::chrono::steady_clock::now();
            projectReducedModel(portVectors);
            romProjectionSec_ += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
            const auto next = nextExpansion(model_, sweepFrequencies, expansionFrequencies);
            maxRelativeResidual_ = next.maxResidual;
            adaptiveConverged_ = maxRelativeResidual_ < residualTolerance;
            if (!adaptiveConverged_ && point + 1 < static_cast<std::size_t>(options_.expansionPointCount)
                && next.frequency > 0.0) {
                expansionFrequencies.push_back(next.frequency);
            }
        }
    }

    romDim_ = static_cast<int>(basis_.size());
    if (romDim_ == 0) {
        throw std::runtime_error("MgaweSweep: global basis collapsed to zero dimension");
    }

    const auto projectionStarted = std::chrono::steady_clock::now();
    if (!options_.adaptive || sweepFrequencies.empty()) projectReducedModel(portVectors);
    basis_.clear();
    romProjectionSec_ += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - projectionStarted).count();
    ready_ = true;
    return romDim_;
}

bool MgaweSweep::appendIfIndependent(std::vector<Complex>&& w) {
    const double pre = fastsweep::norm2(w);
    if (pre < options_.dropTolerance) {
        ++deflatedColumns_;
        return false;
    }
    const double post = fastsweep::modifiedGramSchmidt(w, basis_);
    if (post < options_.dropTolerance * std::max(pre, 1.0)) {
        ++deflatedColumns_;
        return false;
    }
    for (auto& z : w) {
        z /= post;
    }
    basis_.push_back(std::move(w));
    ++retainedColumns_;
    return true;
}

void MgaweSweep::projectReducedModel(const std::vector<std::vector<Complex>>& portVectors) {
    model_.build(project_, affine_, portModeSolver_, basis_, portVectors);
    model_.setPolynomials(polynomials_, options_.useLastPolynomial);
}

std::vector<MgaweSweep::Complex> MgaweSweep::solveReduced(double frequencyHz) const {
    return model_.solveReduced(frequencyHz);
}

SParameterPoint MgaweSweep::evaluate(double frequencyHz) const {
    return model_.evaluate(frequencyHz);
}

std::vector<MgaweSweep::Complex> MgaweSweep::reconstructField(double frequencyHz) const {
    return model_.reconstructField(frequencyHz);
}

const PortMode& MgaweSweep::virtualPortMode(int virtualPortIndex) const {
    return fastsweep::virtualPortMode(portModeSolver_, affine_, virtualPortIndex);
}

}  // namespace fem::sweep
