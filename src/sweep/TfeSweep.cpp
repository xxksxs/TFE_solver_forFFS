#include "bpfem/sweep/TfeSweep.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PolynomialMomentRecurrence.hpp"
#include "bpfem/fastsweep/TfePolynomialModel.hpp"
#include "bpfem/linalg/FactorizedSolveSession.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>

namespace fem::sweep {
namespace {
using Clock = std::chrono::steady_clock;
using tfe::ModeStatus;
using tfe::TfeSystem;
using Vector = TfeSystem::Vector;
using Solution = std::vector<Vector>;
using fastsweep::TfePolynomialModel;

struct Incidence {
    std::vector<std::size_t> input;
    std::vector<Vector> columns;
};

Incidence makeIncidence(const TfeSystem& system, double f, const ProjectDefinition& project,
                        bool all, int requestedPort, int requestedMode) {
    Incidence result;
    const auto& channels = system.channels();
    if (!all && requestedPort < 0) {
        std::set<std::pair<int, int>> seen;
        for (const auto& source : project.modalSources) {
            if (!seen.emplace(source.portId, source.mode).second) {
                throw std::runtime_error("Duplicate AEDT TFE source channel");
            }
            if (source.magnitudeW > 0.0 && std::none_of(channels.begin(), channels.end(), [&](const auto& c) {
                return c.portId == source.portId && c.mode == source.mode;
            })) throw std::runtime_error("Active AEDT source mode is outside the retained TFE spectrum");
        }
    }
    for (std::size_t i = 0; i < channels.size(); ++i) {
        const auto& channel = channels[i];
        if (system.status(i, f) == ModeStatus::Cutoff) {
            // The normalized coordinate is undefined even for a passive mode.
            system.admittance(i, f);
        }
        if (all && system.status(i, f) != ModeStatus::Propagating) continue;
        const auto port = std::find_if(project.ports.begin(), project.ports.end(),
            [&](const auto& p) { return p.id == channel.portId; });
        const auto source = std::find_if(project.modalSources.begin(), project.modalSources.end(),
            [&](const auto& s) { return s.portId == channel.portId && s.mode == channel.mode; });
        const bool selected = requestedPort >= 0
            ? channel.portId == requestedPort && channel.mode == requestedMode
            : !project.modalSources.empty()
                ? source != project.modalSources.end() && source->magnitudeW > 0.0
                : channel.mode == 1 && port != project.ports.end() && port->excited;
        if (!all && !selected) continue;
        if (system.status(i, f) != ModeStatus::Propagating) {
            throw std::runtime_error("TFE requested input mode is below cutoff at " + std::to_string(f) + " Hz");
        }
        Vector incident(channels.size());
        // Integral E x H is one for the textbook normalized pair. Peak
        // phasors give average power |a|^2/2.
        if (all || requestedPort >= 0) {
            incident[i] = 1.0;
        } else {
            const double power = source == project.modalSources.end() ? port->magnitudeW : source->magnitudeW;
            const double phase = source == project.modalSources.end() ? port->phaseDeg : source->phaseDeg;
            incident[i] = std::polar(std::sqrt(2.0 * power), phase * pi / 180.0);
        }
        if (incident[i] == 0.0) continue;
        result.input.push_back(i);
        result.columns.push_back(std::move(incident));
    }
    if (result.columns.empty()) throw std::runtime_error("TFE has no active incidence at " + std::to_string(f) + " Hz");
    return result;
}

struct ReducedBand {
    std::vector<Vector> basis;
    std::array<Vector, 3> matrices;
    std::vector<Vector> unitRhs;

    Solution solve(const TfePolynomialModel::Weights& weights, double frequencyHz) const {
        const auto q = basis.size();
        Vector matrix(q * q);
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t i = 0; i < matrix.size(); ++i) matrix[i] += weights[j] * matrices[j][i];
        }
        Solution states;
        states.reserve(unitRhs.size());
        for (const auto& rhs : unitRhs) {
            auto copy = matrix;
            Vector g = rhs;
            for (auto& coefficient : g) coefficient *= frequencyHz;
            if (!fastsweep::denseSolve(copy, g, static_cast<int>(q))) {
                throw std::runtime_error("Singular TFE-MGAWE reduced system");
            }
            Vector full(basis.front().size());
            for (std::size_t k = 0; k < q; ++k) {
                for (std::size_t i = 0; i < full.size(); ++i) full[i] += g[k] * basis[k][i];
            }
            states.push_back(std::move(full));
        }
        return states;
    }
};

ReducedBand reduce(const TfeSystem& system, const TfePolynomialModel& polynomial,
                   const Incidence& incidence, const MgaweOptions& options,
                   double left, double right, const SweepContext& ctx) {
    ReducedBand model;
    std::vector<std::pair<double, int>> expansions;
    if (options.expansions.empty()) {
        for (int k = 0; k < options.expansionPointCount; ++k) {
            const double ratio = (k + 0.5) / options.expansionPointCount;
            expansions.emplace_back(left + (right - left) * ratio, options.localOrder);
        }
    } else {
        for (const auto& expansion : options.expansions) {
            if (expansion.first >= left && expansion.first <= right) expansions.push_back(expansion);
        }
    }
    if (expansions.empty()) throw std::runtime_error("TFE-MGAWE subband has no expansion point");
    const linalg::SolverConfig config{ctx.linearMaxIterations, ctx.linearTolerance};
    for (const auto& [frequencyHz, order] : expansions) {
        const auto t = polynomial.coordinate(frequencyHz);
        auto a0 = polynomial.matrix(polynomial.weights(t));
        const auto a1 = polynomial.weights(t, 1);
        const auto a2 = polynomial.weights(t, 2);
        linalg::FactorizedSolveSession session(ctx.solver, a0, config);
        const std::vector<fastsweep::PolynomialMomentRecurrence::LinearOperator> operators{
            [&](const auto& x) { return polynomial.apply(a1, x); },
            [&](const auto& x) { return polynomial.apply(a2, x); }};
        for (const auto& input : incidence.columns) {
            // b(f)=4*pi*j*mu0*f*J_alpha*a is exactly linear for fixed incident amplitudes.
            const auto b0 = system.rhs(frequencyHz, input);
            const double width = (right - left) * 0.5;
            auto moments = fastsweep::PolynomialMomentRecurrence::generatePolynomialMoments(
                static_cast<std::size_t>(order),
                [&](std::size_t r) {
                    if (r == 0) return b0;
                    Vector b(b0.size());
                    if (r == 1) {
                        for (std::size_t i = system.interiorCount(); i < b.size(); ++i) {
                            b[i] = b0[i] * (width / frequencyHz);
                        }
                    }
                    return b;
                }, operators, [&](const auto& b) {
                    auto result = session.solve(b);
                    if (result.status != SolveStatus::Converged && ctx.onFrequencyStat) {
                        ctx.onFrequencyStat(frequencyHz, 0.0, 0.0, result);
                    }
                    linalg::requireConverged(result, ctx.solver.name());
                    return std::move(result.field);
                });
            for (auto& moment : moments) {
                const double before = fastsweep::norm2(moment);
                const double after = fastsweep::modifiedGramSchmidt(moment, model.basis);
                if (!(before > 0.0) || after <= options.dropTolerance * before) continue;
                for (auto& value : moment) value /= after;
                model.basis.push_back(std::move(moment));
            }
        }
        ctx.log.info("TFE-MGAWE expansion " + std::to_string(frequencyHz) + " Hz, moments "
                     + std::to_string(order) + ", accumulated rank " + std::to_string(model.basis.size()));
    }
    const auto q = model.basis.size();
    if (q == 0) throw std::runtime_error("Empty TFE-MGAWE basis");
    for (std::size_t k = 0; k < 3; ++k) {
        auto& matrix = model.matrices[k];
        matrix.resize(q * q);
        for (std::size_t j = 0; j < q; ++j) {
            const auto av = polynomial.samples()[k].multiply(model.basis[j]);
            for (std::size_t i = 0; i <= j; ++i) {
                const auto entry = fastsweep::bilinear(model.basis[i], av);
                matrix[i * q + j] = matrix[j * q + i] = entry;
            }
        }
    }
    for (const auto& incident : incidence.columns) {
        Vector rhs(q);
        for (std::size_t j = 0; j < q; ++j) {
            for (std::size_t i = 0; i < incident.size(); ++i) {
                rhs[j] += model.basis[j][system.interiorCount() + i] * incident[i];
            }
            rhs[j] *= std::complex<double>(0.0, 4.0 * pi * mu0);
        }
        model.unitRhs.push_back(std::move(rhs));
    }
    return model;
}

bool sameCutoffRegion(const TfeSystem& system, double first, double next) {
    for (std::size_t i = 0; i < system.channels().size(); ++i) {
        if (system.status(i, first) != system.status(i, next)) return false;
    }
    return true;
}

void appendResults(SweepResult& out, const TfeSystem& system, const Incidence& incidence,
                   const Solution& states, double f, const SweepContext& ctx, bool last) {
    const auto& channels = system.channels();
    SParameterPoint legacy{};
    legacy.frequencyHz = f;
    double power = 0.0;
    for (std::size_t col = 0; col < states.size(); ++col) {
        const auto k = incidence.input[col];
        const auto wave = system.outgoing(states[col], incidence.columns[col]);
        for (std::size_t i = 0; i < channels.size(); ++i) {
            if (system.status(i, f) != ModeStatus::Propagating) continue;
            const auto s = wave[i] / incidence.columns[col][k];
            out.modalScattering.push_back({f, channels[k].portId, channels[k].mode,
                                           channels[i].portId, channels[i].mode, s});
            if (col == 0) {
                power += std::norm(s);
                if (i == k) legacy.s11 = s;
                if (channels[i].portId != channels[k].portId && channels[i].mode == 1) {
                    // The legacy CSV uses the first other physical port.
                    const auto other = std::find_if(channels.begin(), channels.end(),
                        [&](const auto& channel) { return channel.portId != channels[k].portId; });
                    if (i == static_cast<std::size_t>(other - channels.begin())) legacy.s21 = s;
                }
            }
        }
    }
    ctx.log.info("TFE all-propagating output power for first input = " + std::to_string(power)
                 + " at " + std::to_string(f) + " Hz");
    out.points.push_back(legacy);
    if (last || ctx.onFieldSolved) {
        auto edge = system.restoreEdgeField(f, states.front());
        if (ctx.onFieldSolved) ctx.onFieldSolved(f, edge);
        if (last) {
            out.lastEdgeDofs = std::move(edge);
            out.lastFrequencyHz = f;
        }
    }
}
}  // namespace

TfeSweep::TfeSweep(const TfeSystem& system, bool allSources, int sourcePort, int sourceMode)
    : system_(system), allSources_(allSources), sourcePort_(sourcePort), sourceMode_(sourceMode) {}

TfeSweep::TfeSweep(const TfeSystem& system, bool allSources, int sourcePort, int sourceMode,
                   MgaweOptions options, bool proxyDirect, double maxSubbandWidthHz,
                   std::vector<std::pair<double, double>> exactBandsHz)
    : system_(system), allSources_(allSources), sourcePort_(sourcePort), sourceMode_(sourceMode),
      mgawe_(true), proxyDirect_(proxyDirect), maxSubbandWidthHz_(maxSubbandWidthHz),
      exactBandsHz_(std::move(exactBandsHz)),
      options_(std::move(options)) {}

SweepResult TfeSweep::run(const std::vector<double>& frequencies, const SweepContext& ctx) {
    SweepResult out;
    if (mgawe_ && options_.expansions.empty() &&
        (options_.expansionPointCount < 1 || options_.localOrder < 1)) {
        throw std::invalid_argument("TFE-MGAWE requires positive expansion count and order");
    }
    for (const auto& [f, order] : options_.expansions) {
        if (frequencies.empty() || f < frequencies.front() || f > frequencies.back() || order < 1) {
            throw std::invalid_argument("TFE-MGAWE expansion is outside the sweep or has invalid order");
        }
        for (std::size_t i = 0; i < system_.channels().size(); ++i) {
            if (system_.status(i, f) == ModeStatus::Cutoff) {
                throw std::invalid_argument("TFE-MGAWE expansion is on a retained-mode cutoff");
            }
        }
    }
    ctx.solver.rememberPatternForReuse(true);
    const linalg::SolverConfig config{ctx.linearMaxIterations, ctx.linearTolerance};
    const auto useExact = [&](double f) {
        return std::any_of(exactBandsHz_.begin(), exactBandsHz_.end(), [f](const auto& band) {
            return f >= band.first && f <= band.second;
        });
    };
    for (std::size_t begin = 0; begin < frequencies.size();) {
        std::size_t end = begin + 1;
        const auto first = frequencies[begin];
        const auto incidence = makeIncidence(system_, first, ctx.project, allSources_, sourcePort_, sourceMode_);
        while (mgawe_ && end < frequencies.size() &&
               (maxSubbandWidthHz_ <= 0.0 || frequencies[end] - first <= maxSubbandWidthHz_) &&
               sameCutoffRegion(system_, first, frequencies[end])) {
            ++end;
        }
        std::optional<TfePolynomialModel> polynomial;
        ReducedBand model;
        const bool needsProxy = mgawe_ && first < frequencies[end - 1] &&
            std::any_of(frequencies.begin() + begin, frequencies.begin() + end,
                        [&](double f) { return !useExact(f); });
        if (needsProxy) {
            const auto before = Clock::now();
            polynomial.emplace(system_, first, frequencies[end - 1]);
            std::ostringstream message;
            message << std::setprecision(12) << "TFE-MGAWE subband " << first << " to "
                    << frequencies[end - 1] << " Hz; Chebyshev nodes:";
            for (const auto f : polynomial->frequencies()) message << ' ' << f;
            ctx.log.info(message.str());
            if (!proxyDirect_) {
                model = reduce(system_, *polynomial, incidence, options_, first, frequencies[end - 1], ctx);
            }
            out.offlineBuildSec += std::chrono::duration<double>(Clock::now() - before).count();
        } else if (mgawe_ && first == frequencies[end - 1] && !useExact(first)) {
            ctx.log.info("TFE-MGAWE single-point subband; using exact Direct at " + std::to_string(first) + " Hz");
        }
        for (std::size_t index = begin; index < end; ++index) {
            const double f = frequencies[index];
            const bool exactBand = useExact(f);
            const bool useProxy = polynomial.has_value() && !exactBand;
            if (exactBand) {
                ctx.log.info("TFE exact-band Direct at " + std::to_string(f) + " Hz");
            }
            auto before = Clock::now();
            const auto exact = system_.assemble(f);
            std::vector<Vector> rhs;
            for (const auto& input : incidence.columns) rhs.push_back(system_.rhs(f, input));
            const double assembleSec = std::chrono::duration<double>(Clock::now() - before).count();
            before = Clock::now();
            Solution states;
            SolveResult summary;
            summary.status = SolveStatus::Converged;
            const auto collect = [&](std::vector<SolveResult> results) {
                for (auto& result : results) {
                    if (result.status != SolveStatus::Converged && ctx.onFrequencyStat) {
                        ctx.onFrequencyStat(f, assembleSec,
                            std::chrono::duration<double>(Clock::now() - before).count(), result);
                    }
                    linalg::requireConverged(result, ctx.solver.name());
                    summary.iterations = std::max(summary.iterations, result.iterations);
                    summary.residual = std::max(summary.residual, result.residual);
                    summary.absoluteResidual = std::max(summary.absoluteResidual, result.absoluteResidual);
                    states.push_back(std::move(result.field));
                }
            };
            double surrogateResidual = 0.0;
            if (useProxy) {
                const auto weights = polynomial->weights(polynomial->coordinate(f));
                if (proxyDirect_) {
                    const auto proxy = polynomial->matrix(weights);
                    linalg::FactorizedSolveSession session(ctx.solver, proxy, config);
                    collect(session.solveBatch(rhs));
                    surrogateResidual = summary.residual;
                } else {
                    states = model.solve(weights, f);
                }
                if (!proxyDirect_) {
                    for (std::size_t col = 0; col < states.size(); ++col) {
                        auto r = polynomial->apply(weights, states[col]);
                        for (std::size_t i = 0; i < r.size(); ++i) r[i] -= rhs[col][i];
                        const auto residual = linalg::vectorNorm(r) /
                            std::max(linalg::vectorNorm(rhs[col]), 1.0e-300);
                        if (!std::isfinite(residual)) {
                            throw std::runtime_error("Nonfinite TFE quadratic proxy residual");
                        }
                        surrogateResidual = std::max(surrogateResidual, residual);
                    }
                }
            } else {
                linalg::FactorizedSolveSession session(ctx.solver, exact, config);
                collect(session.solveBatch(rhs));
            }
            const double solveSec = std::chrono::duration<double>(Clock::now() - before).count();
            if (useProxy) {
                out.onlineSweepSec += assembleSec + solveSec;
                // The approximation error in the original equation is not a
                // sparse linear solver failure. Report it independently.
                summary.status = SolveStatus::Unchecked;
                summary.residual = 0.0;
                summary.absoluteResidual = 0.0;
                for (std::size_t col = 0; col < states.size(); ++col) {
                    SolveResult checked;
                    checked.field = states[col];
                    linalg::checkResidual(exact, rhs[col], config.tolerance, checked);
                    if (checked.status == SolveStatus::NonFinite) {
                        if (ctx.onFrequencyStat) ctx.onFrequencyStat(f, assembleSec, solveSec, checked);
                        linalg::requireConverged(checked, name());
                    }
                    summary.residual = std::max(summary.residual, checked.residual);
                    summary.absoluteResidual = std::max(summary.absoluteResidual, checked.absoluteResidual);
                }
                std::ostringstream report;
                report << name() << " residual at " << f << " Hz: exact=" << summary.residual
                       << ", polynomial=" << surrogateResidual;
                ctx.log.info(report.str());
            }
            if (ctx.onFrequencyStat) ctx.onFrequencyStat(f, assembleSec, solveSec, summary);
            appendResults(out, system_, incidence, states, f, ctx, index + 1 == frequencies.size());
        }
        begin = end;
    }
    return out;
}

}  // namespace fem::sweep
