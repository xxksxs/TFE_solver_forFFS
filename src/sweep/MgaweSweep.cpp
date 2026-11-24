#include "bpfem/sweep/MgaweSweep.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PolynomialMomentRecurrence.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"
#include "bpfem/linalg/FactorizedSolveSession.hpp"

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
ExpansionChoice nextExpansion(const std::vector<double>& frequencies,
                              const std::vector<double>& expansions,
                              const std::vector<double>& residuals,
                              std::vector<bool>& converged) {
    const int count = static_cast<int>(frequencies.size());
    const double maxResidual = *std::max_element(residuals.begin(), residuals.end());
    if (maxResidual < residualTolerance) return {0.0, maxResidual};

    // Recheck points already covered because the final global ROM is used at every frequency.
    // The paper permits retaining outputs computed with an earlier basis at accepted points.
    for (int i = 0; i < count; ++i) {
        converged[i] = converged[i] && residuals[i] < residualTolerance;
    }
    const auto latest = std::find(frequencies.begin(), frequencies.end(), expansions.back());
    const int center = static_cast<int>(latest - frequencies.begin());
    for (int i = center; i < count && residuals[i] < residualTolerance; ++i) converged[i] = true;
    for (int i = center - 1; i >= 0 && residuals[i] < residualTolerance; --i) converged[i] = true;

    auto used = [&](int i) {
        return std::find(expansions.begin(), expansions.end(), frequencies[i]) != expansions.end();
    };
    int peak = -1;
    for (int i = 0; i < count; ++i) {
        if (!used(i) && !converged[i]
            && (peak < 0 || residuals[i] > residuals[peak])) peak = i;
    }
    if (peak < 0) return {0.0, maxResidual};
    int left = peak;
    int right = peak;
    while (left > 0 && !converged[left - 1] && !used(left - 1)) --left;
    while (right + 1 < count && !converged[right + 1] && !used(right + 1)) ++right;

    const bool touchesLeft = left == 0;
    const bool touchesRight = right == count - 1;
    // The subband endpoints in (19) are grid points, including the band extrema.
    const int leftBoundary = touchesLeft ? 0 : left - 1;
    const int rightBoundary = touchesRight ? count - 1 : right + 1;
    const double leftWeight = touchesLeft ? 0.75 : (touchesRight ? 0.25 : 0.5);
    const int next = std::clamp(static_cast<int>(std::ceil(
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
    if (options_.adaptiveTrueResidual) {
        options_.adaptive = true;
        options_.polynomialProjection = true;
        if (!std::isfinite(options_.trueResidualTolerance) || options_.trueResidualTolerance <= 0.0) {
            throw std::invalid_argument("MGAWE true residual tolerance must be finite and positive");
        }
    }
    if (options_.adaptive && !options_.expansions.empty()) {
        throw std::invalid_argument("MGAWE: explicit expansions are only supported in fixed-point modes");
    }
    if (options_.adaptive && options_.useLastPolynomial) {
        throw std::invalid_argument("MGAWE: --mgawe-2001 and --mgawe-adaptive are mutually exclusive");
    }
    options_.expansionPointCount = std::max(1, options_.expansionPointCount);
    options_.localOrder = std::max(1, options_.localOrder);
    if (!std::isfinite(options_.dropTolerance) || options_.dropTolerance <= 0.0) {
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
        // Equation (18) uses ceil((1 + N) / 2) with one-based indices.
        expansionFrequencies = {frequencies[frequencies.size() / 2]};
    }
    const int dim = build(expansionFrequencies, frequencies, ctx.solver, solverCfg);
    if (options_.adaptive) {
        std::ostringstream summary;
        summary << "MGAWE adaptive stop: "
                << (options_.adaptiveTrueResidual ? engineeringStopReason_
                    : (adaptiveConverged_ ? "residual tolerance" : "configured expansion-point limit"))
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
        if (options_.adaptiveTrueResidual) {
            diag.model = "exact_port_frequency_law_taylor_basis";
            diag.stoppingReason = engineeringStopReason_;
            diag.trueResidualTolerance = options_.trueResidualTolerance;
            diag.residualSamples = std::move(residualSamples_);
            diag.unmetFrequenciesHz = std::move(unmetFrequenciesHz_);
            if (const auto* factorized = dynamic_cast<const linalg::IFactorizedSparseSolver*>(&ctx.solver)) {
                const auto stats = factorized->factorizationStatistics();
                diag.factorizedRhsSolveCount = stats.factorizedRhsSolveCount;
                diag.factorizedSolveCallCount = stats.factorizedSolveCallCount;
                diag.batchRhsMax = stats.batchRhsMax;
            }
        }
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
    engineeringStopReason_.clear();
    residualSamples_.clear();
    unmetFrequenciesHz_.clear();
    if (project_.ports.size() < 2) {
        throw std::runtime_error("MgaweSweep: at least two project ports are required");
    }
    if (expansionFrequencies.empty()) {
        throw std::runtime_error("MgaweSweep: at least one expansion point is required");
    }

    affine_ = assembler_.buildAffineSystem();
    if (affine_.portCoupling.empty()) {
        throw std::runtime_error("MgaweSweep: AffineSystem produced zero port entries");
    }
    if (!affine_.lossless) {
        throw std::runtime_error("MgaweSweep: lossy materials (sigma != 0) are not supported. Use --sweep direct.");
    }

    std::vector<bool> expandable(sweepFrequencies.size(), false);
    if (options_.adaptiveTrueResidual && !sweepFrequencies.empty()) {
        const double center = expansionFrequencies.front();
        std::optional<double> nearest;
        bool anyExcitation = false;
        for (std::size_t i = 0; i < sweepFrequencies.size(); ++i) {
            const double frequency = sweepFrequencies[i];
            const bool excited = fastsweep::hasLosslessExcitation(project_, portModeSolver_, affine_, frequency);
            anyExcitation = anyExcitation || excited;
            expandable[i] = excited && !fastsweep::nearPortCutoff(affine_, frequency);
            if (expandable[i] && (!nearest || std::abs(frequency - center) < std::abs(*nearest - center))) {
                nearest = frequency;
            }
        }
        if (!nearest) {
            if (!anyExcitation) {
                throw std::runtime_error("MGAWE engineering mode: zero excitation RHS at every requested frequency");
            }
            throw std::runtime_error(
                "MGAWE engineering mode: all frequencies with nonzero excitation are inside the port cutoff guard");
        }
        expansionFrequencies.front() = *nearest;
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
    std::vector<bool> converged(sweepFrequencies.size(), false);

    for (std::size_t point = 0; point < expansionFrequencies.size(); ++point) {
        const double expansionHz = expansionFrequencies[point];
        const bool adaptive = options_.adaptive && !sweepFrequencies.empty();
        using Builder = fastsweep::PolynomialPortMomentBuilder;
        const auto variable = options_.polynomialProjection ? Builder::Variable::WavenumberShift
                                                            : Builder::Variable::LambdaShift;
        const int requestedOrder = !options_.expansions.empty()
            ? options_.expansions.at(point).second : options_.localOrder;
        // Keep at least three Taylor coefficients for the exact quadratic mass term online.
        // Early exhaustion limits the number of moments, independently of this truncation.
        auto series = Builder::buildTaylorSeries(project_, portModeSolver_, affine_, expansionHz,
            options_.polynomialProjection ? std::max(3, requestedOrder) : requestedOrder, variable);
        std::size_t solvedMoments = 0;
        {
            Vector rhs0;
            const auto a0 = assembler_.assemble(expansionHz, rhs0);
            using Recurrence = fastsweep::PolynomialMomentRecurrence;
            std::vector<Recurrence::LinearOperator> operators;
            for (int r = 1; r < requestedOrder; ++r) {
                operators.emplace_back([&, r](const Vector& x) {
                    return Builder::applyMatrixCoefficient(affine_, series, portVectors,
                        static_cast<std::size_t>(r), x);
                });
            }
            linalg::FactorizedSolveSession session(solver, a0, solverConfig);
            Recurrence::generatePolynomialMoments(static_cast<std::size_t>(requestedOrder),
                [&](std::size_t r) {
                    return r == 0 ? rhs0 : Builder::rhsCoefficient(series, portVectors, rhs0.size(), r);
                },
                operators, [&](const Vector& rhs) {
                    auto result = session.solve(rhs);
                    if (result.status != SolveStatus::Converged) {
                        const auto context = std::string(solver.name()) + " at expansion "
                            + std::to_string(expansionHz) + " Hz, moment " + std::to_string(solvedMoments);
                        linalg::requireConverged(result, context.c_str());
                    }
                    ++solvedMoments;
                    const auto start = std::chrono::steady_clock::now();
                    appendIfIndependent(Vector(result.field));
                    orthogonalizationSec_ += std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - start).count();
                    return std::move(result.field);
                }, [&](const Vector& rhs) {
                    return !adaptive || options_.adaptiveTrueResidual || basis_.empty()
                        || !localSpaceExhausted(a0, basis_, rhs);
                });
        }
        localOrders_.push_back(static_cast<int>(solvedMoments));
        if (options_.polynomialProjection) {
            polynomials_.push_back(std::move(series));
        }
        if (adaptive && !basis_.empty()) {
            const auto start = std::chrono::steady_clock::now();
            projectReducedModel(portVectors);
            romProjectionSec_ += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
            if (options_.adaptiveTrueResidual) {
                const auto next = assessTrueResiduals(sweepFrequencies, expansionFrequencies, expandable);
                if (next) expansionFrequencies.push_back(*next);
                continue;
            }
            std::vector<double> residuals;
            residuals.reserve(sweepFrequencies.size());
            for (double frequency : sweepFrequencies) {
                const auto* series = polynomialAt(frequency);
                const double residual = options_.polynomialProjection && !series
                    ? std::numeric_limits<double>::infinity() : model_.relativeResidual(frequency, series);
                residuals.push_back(std::isfinite(residual) ? residual : std::numeric_limits<double>::infinity());
            }
            const auto next = nextExpansion(sweepFrequencies, expansionFrequencies, residuals, converged);
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

std::optional<double> MgaweSweep::assessTrueResiduals(const std::vector<double>& frequencies,
    const std::vector<double>& expansions, const std::vector<bool>& expandable) {
    residualSamples_.clear();
    unmetFrequenciesHz_.clear();
    maxRelativeResidual_ = 0.0;
    std::optional<double> next;
    double worstCandidate = -1.0;
    for (std::size_t index = 0; index < frequencies.size(); ++index) {
        const double frequency = frequencies[index];
        const bool hasExcitation = expandable[index]
            || fastsweep::hasLosslessExcitation(project_, portModeSolver_, affine_, frequency);
        const auto field = hasExcitation ? model_.reconstructField(frequency) : Vector(affine_.K.size());
        auto sample = fastsweep::losslessResiduals(project_, portModeSolver_, affine_,
            frequency, field, polynomialAt(frequency));
        if (!sample.zeroRhs) {
            const double relative = sample.trueRelative.value_or(std::numeric_limits<double>::infinity());
            maxRelativeResidual_ = std::max(maxRelativeResidual_, relative);
            if (relative > options_.trueResidualTolerance) {
                unmetFrequenciesHz_.push_back(frequency);
                const bool unused = std::find(expansions.begin(), expansions.end(), frequency) == expansions.end();
                if (expandable[index] && unused && relative > worstCandidate) {
                    next = frequency;
                    worstCandidate = relative;
                }
            }
        }
        residualSamples_.push_back(sample);
    }
    adaptiveConverged_ = unmetFrequenciesHz_.empty();
    if (adaptiveConverged_) {
        engineeringStopReason_ = "true_residual_tolerance";
        return std::nullopt;
    }
    if (expansions.size() >= static_cast<std::size_t>(options_.expansionPointCount)) {
        engineeringStopReason_ = "expansion_budget_exhausted";
        return std::nullopt;
    }
    if (!next) engineeringStopReason_ = "no_eligible_unmet_frequency";
    return next;
}

bool MgaweSweep::appendIfIndependent(std::vector<Complex>&& w) {
    const double pre = fastsweep::norm2(w);
    if (!std::isfinite(pre)) {
        throw std::runtime_error("MGAWE: a local moment has a non-finite norm");
    }
    if (pre == 0.0) {
        ++deflatedColumns_;
        return false;
    }
    // Test the fraction outside the current space. The raw moment norm depends on
    // the excitation amplitude and on the scale of the Taylor expansion variable.
    for (auto& z : w) {
        z /= pre;
    }
    const double post = fastsweep::modifiedGramSchmidt(w, basis_);
    if (post <= options_.dropTolerance) {
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
    using Projection = fastsweep::GalerkinReducedModel::Projection;
    model_.build(project_, affine_, portModeSolver_, basis_, portVectors,
                 options_.useLastPolynomial ? Projection::Hermitian : Projection::Transpose);
}

const MgaweSweep::PortTaylorSeries* MgaweSweep::polynomialAt(double frequencyHz) const {
    const double k = 2.0 * pi * frequencyHz / c0;
    const double lambda = k * k;
    const PortTaylorSeries* nearest = nullptr;
    for (const auto& candidate : polynomials_) {
        if (options_.useLastPolynomial && &candidate != &polynomials_.back()) continue;
        const bool sameBranch = std::all_of(affine_.portCutoffSquared.begin(),
            affine_.portCutoffSquared.end(), [&](double cutoff) {
                return (lambda > cutoff) == (candidate.lambda0 > cutoff);
            });
        if (sameBranch && (!nearest || std::abs(frequencyHz - candidate.expansionHz)
                                       < std::abs(frequencyHz - nearest->expansionHz))) {
            nearest = &candidate;
        }
    }
    return nearest;
}

const MgaweSweep::PortTaylorSeries* MgaweSweep::requirePolynomialAt(double frequencyHz) const {
    if (options_.adaptiveTrueResidual) return nullptr;
    const auto* series = polynomialAt(frequencyHz);
    if (options_.polynomialProjection && !series) {
        throw std::runtime_error("MGAWE: no polynomial expansion on the port cutoff branch at "
                                 + std::to_string(frequencyHz) + " Hz");
    }
    return series;
}

SParameterPoint MgaweSweep::evaluate(double frequencyHz) const {
    if (!ready_) return SParameterPoint{frequencyHz};
    if (options_.adaptiveTrueResidual
        && !fastsweep::hasLosslessExcitation(project_, portModeSolver_, affine_, frequencyHz)) {
        return SParameterPoint{frequencyHz};
    }
    return model_.evaluate(frequencyHz, requirePolynomialAt(frequencyHz));
}

std::vector<MgaweSweep::Complex> MgaweSweep::reconstructField(double frequencyHz) const {
    if (options_.adaptiveTrueResidual
        && !fastsweep::hasLosslessExcitation(project_, portModeSolver_, affine_, frequencyHz)) {
        return Vector(affine_.K.size());
    }
    return model_.reconstructField(frequencyHz, requirePolynomialAt(frequencyHz));
}

}  // namespace fem::sweep
