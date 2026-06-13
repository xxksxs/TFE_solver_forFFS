#include "bpfem/sweep/WcaweSweep.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/PolynomialMomentRecurrence.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fem::sweep {

namespace {

using Complex = std::complex<double>;

std::vector<Complex> liftPortVector(const std::vector<std::pair<int, double>>& sparse,
                                    std::size_t n) {
    std::vector<Complex> out(n, Complex(0.0, 0.0));
    for (const auto& [i, w] : sparse) {
        if (i >= 0 && static_cast<std::size_t>(i) < n) {
            out[static_cast<std::size_t>(i)] = w;
        }
    }
    return out;
}

Complex bilinear(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    Complex s(0.0, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        s += a[i] * b[i];
    }
    return s;
}

Complex hdot(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    Complex s(0.0, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        s += std::conj(a[i]) * b[i];
    }
    return s;
}

double norm2(const std::vector<Complex>& a) {
    double s = 0.0;
    for (const auto& z : a) {
        s += std::norm(z);
    }
    return std::sqrt(s);
}

double binomialCoefficient(double alpha, int order) {
    double out = 1.0;
    for (int j = 1; j <= order; ++j) {
        out *= (alpha - static_cast<double>(j - 1)) / static_cast<double>(j);
    }
    return out;
}

std::vector<double> shiftedPowerSeries(double base, double alpha, int count) {
    if (base <= 0.0) {
        throw std::runtime_error("WcaweSweep: expansion series base must be positive");
    }
    std::vector<double> coeffs(static_cast<std::size_t>(count), 0.0);
    double invPower = 1.0;
    for (int r = 0; r < count; ++r) {
        coeffs[static_cast<std::size_t>(r)] = binomialCoefficient(alpha, r) * invPower;
        invPower /= base;
    }
    return coeffs;
}

std::vector<double> convolve(const std::vector<double>& a,
                             const std::vector<double>& b,
                             int count) {
    std::vector<double> out(static_cast<std::size_t>(count), 0.0);
    for (int i = 0; i < count; ++i) {
        double sum = 0.0;
        for (int j = 0; j <= i; ++j) {
            sum += a[static_cast<std::size_t>(j)] * b[static_cast<std::size_t>(i - j)];
        }
        out[static_cast<std::size_t>(i)] = sum;
    }
    return out;
}

double powInt(double x, int n) {
    double out = 1.0;
    for (int i = 0; i < n; ++i) {
        out *= x;
    }
    return out;
}

bool denseSolve(std::vector<Complex>& A, std::vector<Complex>& rhs, int n) {
    auto idx = [n](int r, int c) {
        return static_cast<std::size_t>(r) * static_cast<std::size_t>(n)
            + static_cast<std::size_t>(c);
    };

    for (int k = 0; k < n; ++k) {
        int pivotRow = k;
        double pivotMag = std::abs(A[idx(k, k)]);
        for (int r = k + 1; r < n; ++r) {
            const double m = std::abs(A[idx(r, k)]);
            if (m > pivotMag) {
                pivotMag = m;
                pivotRow = r;
            }
        }
        if (pivotMag < 1.0e-30) {
            return false;
        }
        if (pivotRow != k) {
            for (int c = 0; c < n; ++c) {
                std::swap(A[idx(k, c)], A[idx(pivotRow, c)]);
            }
            std::swap(rhs[static_cast<std::size_t>(k)],
                      rhs[static_cast<std::size_t>(pivotRow)]);
        }

        const Complex pivot = A[idx(k, k)];
        for (int r = k + 1; r < n; ++r) {
            const Complex factor = A[idx(r, k)] / pivot;
            A[idx(r, k)] = factor;
            for (int c = k + 1; c < n; ++c) {
                A[idx(r, c)] -= factor * A[idx(k, c)];
            }
            rhs[static_cast<std::size_t>(r)] -= factor * rhs[static_cast<std::size_t>(k)];
        }
    }

    for (int r = n - 1; r >= 0; --r) {
        Complex sum = rhs[static_cast<std::size_t>(r)];
        for (int c = r + 1; c < n; ++c) {
            sum -= A[idx(r, c)] * rhs[static_cast<std::size_t>(c)];
        }
        rhs[static_cast<std::size_t>(r)] = sum / A[idx(r, r)];
    }
    return true;
}

}  // namespace

WcaweSweep::WcaweSweep(const ProjectDefinition& project,
                       const FEMAssembler& assembler,
                       const PortModeSolver& portModeSolver,
                       WcaweOptions options)
    : project_(project),
      assembler_(assembler),
      portModeSolver_(portModeSolver),
      options_(options) {
    options_.order = std::max(1, options_.order);
    if (options_.dropTolerance <= 0.0) {
        options_.dropTolerance = 1.0e-12;
    }
}

SweepResult WcaweSweep::run(const std::vector<double>& frequencies, const SweepContext& ctx) {
    const double expansionHz = options_.expansionFrequencyHz > 0.0
        ? options_.expansionFrequencyHz
        : (frequencies.empty()
            ? 0.5 * (project_.sweep.startHz + project_.sweep.endHz)
            : 0.5 * (frequencies.front() + frequencies.back()));

    ctx.log.info("WCAWE expansion frequency: " + std::to_string(expansionHz / 1.0e9) + " GHz");
    ctx.log.info("WCAWE target basis order: " + std::to_string(options_.order));
    std::ostringstream dropTol;
    dropTol << std::scientific << std::setprecision(3) << options_.dropTolerance;
    ctx.log.info("WCAWE MGS/R diagonal drop tolerance: " + dropTol.str());

    linalg::SolverConfig solverCfg;
    solverCfg.maxIterations = ctx.linearMaxIterations;
    solverCfg.tolerance = ctx.linearTolerance;
    const int dim = buildOffline(expansionHz, ctx.solver, solverCfg);
    ctx.log.info("WCAWE ROM dimension: " + std::to_string(dim));
    ctx.log.info("WCAWE deflated candidate columns: " + std::to_string(deflatedColumns_));

    if (!ctx.outputDirectory.empty()) {
        const auto conditionPath = ctx.outputDirectory / "basis_condition.csv";
        if (writeBasisConditionCsv(conditionPath)) {
            ctx.log.info("Wrote " + conditionPath.string());
        } else {
            ctx.log.warn("Failed to write " + conditionPath.string());
        }
    }

    SweepResult out;
    out.points.reserve(frequencies.size());
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        const double f = frequencies[i];
        if (i == 0 || (i + 1) % 50 == 0 || i + 1 == frequencies.size()) {
            ctx.log.info("WCAWE evaluating point " + std::to_string(i + 1) + "/"
                         + std::to_string(frequencies.size())
                         + " @ " + std::to_string(f / 1.0e9) + " GHz");
        }
        out.points.push_back(evaluate(f));
        out.lastFrequencyHz = f;
        if (ctx.onFieldSolved) {
            ctx.onFieldSolved(f, reconstructField(f));
        }
    }
    if (!frequencies.empty()) {
        out.lastEdgeDofs = reconstructField(frequencies.back());
    }
    return out;
}

int WcaweSweep::buildOffline(double expansionFrequencyHz,
                             linalg::ISparseSolver& solver,
                             const linalg::SolverConfig& solverConfig) {
    if (project_.ports.size() < 2) {
        throw std::runtime_error("WcaweSweep: at least two project ports are required");
    }

    expansionFrequencyHz_ = expansionFrequencyHz;
    numProjectPorts_ = static_cast<int>(project_.ports.size());
    affine_ = assembler_.buildAffineSystem();
    fullDim_ = static_cast<int>(affine_.K.size());
    numVirtualPorts_ = static_cast<int>(affine_.portCoupling.size());
    if (numVirtualPorts_ == 0) {
        throw std::runtime_error("WcaweSweep: AffineSystem produced zero port entries");
    }
    if (!affine_.lossless) {
        throw std::runtime_error("WcaweSweep: lossy materials (sigma != 0) are not supported. Use --sweep direct.");
    }

    std::vector<std::vector<Complex>> portVectors;
    portVectors.reserve(static_cast<std::size_t>(numVirtualPorts_));
    for (const auto& coupling : affine_.portCoupling) {
        portVectors.push_back(liftPortVector(coupling, static_cast<std::size_t>(fullDim_)));
    }

    basis_.clear();
    triangularR_.clear();
    rDiagonalAbs_.clear();
    conditionRecords_.clear();
    deflatedColumns_ = 0;
    basis_.reserve(static_cast<std::size_t>(options_.order));

    const auto aweMoments = buildAweMoments(expansionFrequencyHz, solver, solverConfig, portVectors);
    for (const auto& moment : aweMoments) {
        appendWellConditionedColumn(moment);
    }

    romDim_ = static_cast<int>(basis_.size());
    if (romDim_ == 0) {
        throw std::runtime_error("WcaweSweep: well-conditioned basis collapsed to zero dimension");
    }

    projectReducedModel(portVectors);
    ready_ = true;
    return romDim_;
}

std::vector<std::vector<WcaweSweep::Complex>> WcaweSweep::buildAweMoments(
    double expansionFrequencyHz,
    linalg::ISparseSolver& solver,
    const linalg::SolverConfig& solverConfig,
    const std::vector<std::vector<Complex>>& portVectors) const {
    const double kExpansion = 2.0 * pi * expansionFrequencyHz / c0;
    const double lambda0 = kExpansion * kExpansion;
    const double lambdaScale = std::max(lambda0, 1.0);

    std::vector<Complex> rhs0;
    SparseMatrix A0 = assembler_.assemble(expansionFrequencyHz, rhs0);
    if (rhs0.empty()) {
        throw std::runtime_error("WcaweSweep: expansion RHS is empty");
    }

    auto solveAtExpansion = [&](const std::vector<Complex>& rhs) {
        SolveResult result = solver.solve(A0, rhs, solverConfig);
        return result.field;
    };

    const int momentCount = std::max(1, options_.order);
    const std::size_t fullDim = rhs0.size();
    std::vector<std::vector<double>> admittanceCoeffs(
        static_cast<std::size_t>(numVirtualPorts_),
        std::vector<double>(static_cast<std::size_t>(momentCount), 0.0));
    std::vector<std::vector<Complex>> rhsCoefficients(
        static_cast<std::size_t>(momentCount),
        std::vector<Complex>(fullDim, Complex(0.0, 0.0)));

    for (int v = 0; v < numVirtualPorts_; ++v) {
        const double kc2 = affine_.portCutoffSquared[static_cast<std::size_t>(v)];
        const double betaBase = lambda0 - kc2;
        const bool propagatingAtExpansion = betaBase > 0.0;
        const double seriesBase = propagatingAtExpansion ? betaBase : lambda0;
        const double scalar0 = std::sqrt(seriesBase);
        const auto sqrtSeries = shiftedPowerSeries(seriesBase, 0.5, momentCount);
        for (int r = 0; r < momentCount; ++r) {
            admittanceCoeffs[static_cast<std::size_t>(v)][static_cast<std::size_t>(r)] =
                scalar0 * sqrtSeries[static_cast<std::size_t>(r)] * powInt(lambdaScale, r);
        }

        if (!affine_.isExcitationMode[static_cast<std::size_t>(v)]) {
            continue;
        }

        const PortMode& mode = virtualPortMode(v);
        const double s0 = powerNormalizationFactor(mode, expansionFrequencyHz);
        if (s0 <= 0.0) {
            continue;
        }
        const int projectPort = affine_.projectPortIndex[static_cast<std::size_t>(v)];
        const auto& port = project_.ports[static_cast<std::size_t>(projectPort)];
        const Complex incident0 = std::polar(
            std::sqrt(std::max(port.magnitudeW, 0.0)) * s0,
            port.phaseDeg * pi / 180.0);
        const Complex rhsFactor0 = Complex(0.0, 2.0 * admittanceCoeffs[static_cast<std::size_t>(v)][0])
            * incident0;

        const auto betaQuarter = shiftedPowerSeries(seriesBase, 0.25, momentCount);
        const auto lambdaQuarter = shiftedPowerSeries(lambda0, 0.25, momentCount);
        const auto rhsRatio = propagatingAtExpansion
            ? convolve(betaQuarter, lambdaQuarter, momentCount)
            : std::vector<double>(static_cast<std::size_t>(momentCount), 0.0);
        for (int r = 0; r < momentCount; ++r) {
            const Complex coeff = rhsFactor0
                * rhsRatio[static_cast<std::size_t>(r)]
                * powInt(lambdaScale, r);
            auto& br = rhsCoefficients[static_cast<std::size_t>(r)];
            const auto& m = portVectors[static_cast<std::size_t>(v)];
            for (std::size_t i = 0; i < fullDim; ++i) {
                br[i] += coeff * m[i];
            }
        }
    }
    rhsCoefficients.front() = rhs0;

    std::vector<fastsweep::PolynomialMomentRecurrence::LinearOperator> matrixCoefficientOperators;
    matrixCoefficientOperators.reserve(static_cast<std::size_t>(std::max(0, momentCount - 1)));
    for (int r = 1; r < momentCount; ++r) {
        matrixCoefficientOperators.push_back(
            [this, r, lambdaScale, &admittanceCoeffs, &portVectors](const std::vector<Complex>& x) {
                std::vector<Complex> out(x.size(), Complex(0.0, 0.0));
                if (r == 1) {
                    std::vector<Complex> mx = affine_.M.multiply(x);
                    for (std::size_t i = 0; i < out.size(); ++i) {
                        out[i] -= lambdaScale * mx[i];
                    }
                }
                for (int v = 0; v < numVirtualPorts_; ++v) {
                    const double coeff =
                        admittanceCoeffs[static_cast<std::size_t>(v)][static_cast<std::size_t>(r)];
                    if (coeff == 0.0) {
                        continue;
                    }
                    const auto& m = portVectors[static_cast<std::size_t>(v)];
                    const Complex projection = bilinear(m, x);
                    const Complex scale(0.0, coeff);
                    for (std::size_t i = 0; i < out.size(); ++i) {
                        out[i] += scale * projection * m[i];
                    }
                }
                return out;
            });
    }

    return fastsweep::PolynomialMomentRecurrence::generatePolynomialMoments(
        rhsCoefficients, matrixCoefficientOperators, solveAtExpansion);
}

bool WcaweSweep::appendWellConditionedColumn(const std::vector<Complex>& moment) {
    std::vector<Complex> w = moment;
    const double pre = norm2(w);
    if (pre < options_.dropTolerance) {
        ++deflatedColumns_;
        return false;
    }

    const int oldDim = static_cast<int>(basis_.size());
    std::vector<Complex> coeffs(static_cast<std::size_t>(oldDim), Complex(0.0, 0.0));
    for (int pass = 0; pass < 2; ++pass) {
        for (int j = 0; j < oldDim; ++j) {
            const Complex h = hdot(basis_[static_cast<std::size_t>(j)], w);
            coeffs[static_cast<std::size_t>(j)] += h;
            for (std::size_t i = 0; i < w.size(); ++i) {
                w[i] -= h * basis_[static_cast<std::size_t>(j)][i];
            }
        }
    }

    const double post = norm2(w);
    if (post < options_.dropTolerance * std::max(pre, 1.0)) {
        ++deflatedColumns_;
        return false;
    }

    for (auto& z : w) {
        z /= post;
    }
    basis_.push_back(std::move(w));
    rDiagonalAbs_.push_back(post);

    const int newDim = oldDim + 1;
    std::vector<Complex> newR(static_cast<std::size_t>(newDim) * static_cast<std::size_t>(newDim),
                              Complex(0.0, 0.0));
    for (int r = 0; r < oldDim; ++r) {
        for (int c = 0; c < oldDim; ++c) {
            newR[static_cast<std::size_t>(r) * static_cast<std::size_t>(newDim)
                 + static_cast<std::size_t>(c)] =
                triangularR_[static_cast<std::size_t>(r) * static_cast<std::size_t>(oldDim)
                             + static_cast<std::size_t>(c)];
        }
    }
    for (int r = 0; r < oldDim; ++r) {
        newR[static_cast<std::size_t>(r) * static_cast<std::size_t>(newDim)
             + static_cast<std::size_t>(oldDim)] = coeffs[static_cast<std::size_t>(r)];
    }
    newR[static_cast<std::size_t>(oldDim) * static_cast<std::size_t>(newDim)
         + static_cast<std::size_t>(oldDim)] = Complex(post, 0.0);
    triangularR_ = std::move(newR);

    double minDiag = rDiagonalAbs_.front();
    double maxDiag = rDiagonalAbs_.front();
    for (double d : rDiagonalAbs_) {
        minDiag = std::min(minDiag, d);
        maxDiag = std::max(maxDiag, d);
    }

    double orthError = 0.0;
    for (int i = 0; i < newDim; ++i) {
        for (int j = 0; j < newDim; ++j) {
            const Complex expected = (i == j) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
            orthError = std::max(orthError, std::abs(hdot(basis_[static_cast<std::size_t>(i)],
                                                          basis_[static_cast<std::size_t>(j)]) - expected));
        }
    }

    ConditionRecord rec;
    rec.order = newDim;
    rec.aweConditionProxy = maxDiag / std::max(minDiag, 1.0e-300);
    rec.wcaweConditionProxy = 1.0 + orthError;
    rec.rDiagonalAbs = post;
    rec.orthogonalityError = orthError;
    conditionRecords_.push_back(rec);
    return true;
}

void WcaweSweep::projectReducedModel(const std::vector<std::vector<Complex>>& portVectors) {
    Ktilde_.assign(static_cast<std::size_t>(romDim_) * static_cast<std::size_t>(romDim_),
                   Complex(0.0, 0.0));
    Mtilde_.assign(static_cast<std::size_t>(romDim_) * static_cast<std::size_t>(romDim_),
                   Complex(0.0, 0.0));

    for (int j = 0; j < romDim_; ++j) {
        const auto Kvj = affine_.K.multiply(basis_[static_cast<std::size_t>(j)]);
        const auto Mvj = affine_.M.multiply(basis_[static_cast<std::size_t>(j)]);
        for (int i = 0; i < romDim_; ++i) {
            Ktilde_[static_cast<std::size_t>(i) * static_cast<std::size_t>(romDim_)
                    + static_cast<std::size_t>(j)] =
                bilinear(basis_[static_cast<std::size_t>(i)], Kvj);
            Mtilde_[static_cast<std::size_t>(i) * static_cast<std::size_t>(romDim_)
                    + static_cast<std::size_t>(j)] =
                bilinear(basis_[static_cast<std::size_t>(i)], Mvj);
        }
    }

    portModeReduced_.assign(static_cast<std::size_t>(numVirtualPorts_), {});
    for (int v = 0; v < numVirtualPorts_; ++v) {
        portModeReduced_[static_cast<std::size_t>(v)].assign(static_cast<std::size_t>(romDim_),
                                                              Complex(0.0, 0.0));
        for (int i = 0; i < romDim_; ++i) {
            portModeReduced_[static_cast<std::size_t>(v)][static_cast<std::size_t>(i)] =
                bilinear(basis_[static_cast<std::size_t>(i)], portVectors[static_cast<std::size_t>(v)]);
        }
    }
}

std::vector<WcaweSweep::Complex> WcaweSweep::solveReduced(double frequencyHz) const {
    if (!ready_) {
        throw std::runtime_error("WcaweSweep::solveReduced called before buildOffline");
    }

    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double k0Sq = k0 * k0;
    std::vector<double> beta(static_cast<std::size_t>(numVirtualPorts_), 0.0);
    std::vector<Complex> incident(static_cast<std::size_t>(numVirtualPorts_), Complex(0.0, 0.0));

    for (int v = 0; v < numVirtualPorts_; ++v) {
        const double kc2 = affine_.portCutoffSquared[static_cast<std::size_t>(v)];
        const double betaSq = k0Sq - kc2;
        beta[static_cast<std::size_t>(v)] = (betaSq > 0.0) ? std::sqrt(betaSq) : 0.0;

        if (affine_.isExcitationMode[static_cast<std::size_t>(v)]) {
            const PortMode& pm = virtualPortMode(v);
            const double s = powerNormalizationFactor(pm, frequencyHz);
            if (s > 0.0) {
                const auto& port = project_.ports[static_cast<std::size_t>(
                    affine_.projectPortIndex[static_cast<std::size_t>(v)])];
                incident[static_cast<std::size_t>(v)] = std::polar(
                    std::sqrt(std::max(port.magnitudeW, 0.0)) * s,
                    port.phaseDeg * pi / 180.0);
            }
        }
    }

    const int q = romDim_;
    std::vector<Complex> Atilde(static_cast<std::size_t>(q) * static_cast<std::size_t>(q),
                                Complex(0.0, 0.0));
    for (std::size_t k = 0; k < Atilde.size(); ++k) {
        Atilde[k] = Ktilde_[k] - k0Sq * Mtilde_[k];
    }
    for (int v = 0; v < numVirtualPorts_; ++v) {
        const double bv = beta[static_cast<std::size_t>(v)];
        if (bv <= 0.0) {
            continue;
        }
        const Complex jbeta(0.0, bv);
        const auto& vp = portModeReduced_[static_cast<std::size_t>(v)];
        for (int i = 0; i < q; ++i) {
            const Complex jbVi = jbeta * vp[static_cast<std::size_t>(i)];
            for (int j = 0; j < q; ++j) {
                Atilde[static_cast<std::size_t>(i) * static_cast<std::size_t>(q)
                       + static_cast<std::size_t>(j)] +=
                    jbVi * vp[static_cast<std::size_t>(j)];
            }
        }
    }

    std::vector<Complex> rhs(static_cast<std::size_t>(q), Complex(0.0, 0.0));
    for (int v = 0; v < numVirtualPorts_; ++v) {
        const double bv = beta[static_cast<std::size_t>(v)];
        if (bv <= 0.0 || incident[static_cast<std::size_t>(v)] == Complex(0.0, 0.0)) {
            continue;
        }
        const Complex factor = Complex(0.0, 2.0) * bv * incident[static_cast<std::size_t>(v)];
        const auto& vp = portModeReduced_[static_cast<std::size_t>(v)];
        for (int i = 0; i < q; ++i) {
            rhs[static_cast<std::size_t>(i)] += factor * vp[static_cast<std::size_t>(i)];
        }
    }

    std::vector<Complex> xTilde = rhs;
    if (!denseSolve(Atilde, xTilde, q)) {
        throw std::runtime_error("WcaweSweep::solveReduced: reduced system singular at "
                                 + std::to_string(frequencyHz) + " Hz");
    }
    return xTilde;
}

SParameterPoint WcaweSweep::evaluate(double frequencyHz) const {
    SParameterPoint sp;
    sp.frequencyHz = frequencyHz;
    if (!ready_ || numProjectPorts_ < 2) {
        return sp;
    }

    const auto xTilde = solveReduced(frequencyHz);
    std::vector<double> sNorm(static_cast<std::size_t>(numVirtualPorts_), 0.0);
    std::vector<int> dominantVirtualByProject(static_cast<std::size_t>(numProjectPorts_), -1);

    for (int v = 0; v < numVirtualPorts_; ++v) {
        const int faceId = affine_.portFaceIds[static_cast<std::size_t>(v)];
        const auto* multi = portModeSolver_.multiMode(faceId);
        if (multi != nullptr) {
            int modeIdx = 0;
            for (int u = 0; u < v; ++u) {
                if (affine_.portFaceIds[static_cast<std::size_t>(u)] == faceId) {
                    ++modeIdx;
                }
            }
            const PortMode& pm = multi->modes[static_cast<std::size_t>(modeIdx)];
            sNorm[static_cast<std::size_t>(v)] = powerNormalizationFactor(pm, frequencyHz);
            if (modeIdx == multi->excitationModeIndex) {
                const int prj = affine_.projectPortIndex[static_cast<std::size_t>(v)];
                dominantVirtualByProject[static_cast<std::size_t>(prj)] = v;
            }
        } else {
            const PortMode& pm = portModeSolver_.solve(faceId);
            sNorm[static_cast<std::size_t>(v)] = powerNormalizationFactor(pm, frequencyHz);
            const int prj = affine_.projectPortIndex[static_cast<std::size_t>(v)];
            dominantVirtualByProject[static_cast<std::size_t>(prj)] = v;
        }
    }

    auto portProjection = [&](int virtualIdx) {
        const auto& vp = portModeReduced_[static_cast<std::size_t>(virtualIdx)];
        Complex sum(0.0, 0.0);
        for (int i = 0; i < romDim_; ++i) {
            sum += vp[static_cast<std::size_t>(i)] * xTilde[static_cast<std::size_t>(i)];
        }
        return sum;
    };

    const int vIn = dominantVirtualByProject[0];
    const int vOut = (numProjectPorts_ >= 2) ? dominantVirtualByProject[1] : -1;
    if (vIn < 0) {
        return sp;
    }
    const double sIn = sNorm[static_cast<std::size_t>(vIn)];
    const double sOut = (vOut >= 0) ? sNorm[static_cast<std::size_t>(vOut)] : 0.0;
    if (sIn <= 0.0) {
        return sp;
    }

    const Complex bIn = portProjection(vIn) / sIn;
    const Complex bOut = (vOut >= 0 && sOut > 0.0)
        ? portProjection(vOut) / sOut
        : Complex(0.0, 0.0);
    const auto& port0 = project_.ports[0];
    const Complex incidentExtract = port0.excited
        ? std::polar(std::sqrt(std::max(port0.magnitudeW, 0.0)),
                     port0.phaseDeg * pi / 180.0)
        : Complex(1.0, 0.0);
    if (incidentExtract == Complex(0.0, 0.0)) {
        return sp;
    }
    sp.s11 = (bIn - incidentExtract) / incidentExtract;
    sp.s21 = bOut / incidentExtract;
    return sp;
}

std::vector<WcaweSweep::Complex> WcaweSweep::reconstructField(double frequencyHz) const {
    const auto xTilde = solveReduced(frequencyHz);
    std::vector<Complex> field(static_cast<std::size_t>(fullDim_), Complex(0.0, 0.0));
    for (int j = 0; j < romDim_; ++j) {
        const auto& v = basis_[static_cast<std::size_t>(j)];
        const Complex coeff = xTilde[static_cast<std::size_t>(j)];
        for (int i = 0; i < fullDim_; ++i) {
            field[static_cast<std::size_t>(i)] += coeff * v[static_cast<std::size_t>(i)];
        }
    }
    return field;
}

bool WcaweSweep::writeBasisConditionCsv(const std::filesystem::path& path) const {
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out.good()) {
        return false;
    }
    out << "order,awe_condition_proxy,wcawe_condition_proxy,r_diagonal_abs,orthogonality_error\n";
    out << std::setprecision(17);
    for (const auto& r : conditionRecords_) {
        out << r.order << ','
            << r.aweConditionProxy << ','
            << r.wcaweConditionProxy << ','
            << r.rDiagonalAbs << ','
            << r.orthogonalityError << '\n';
    }
    out.flush();
    return out.good();
}

const PortMode& WcaweSweep::virtualPortMode(int virtualPortIndex) const {
    const int faceId = affine_.portFaceIds[static_cast<std::size_t>(virtualPortIndex)];
    if (const auto* multi = portModeSolver_.multiMode(faceId); multi != nullptr) {
        int modeIndex = 0;
        for (int i = 0; i < virtualPortIndex; ++i) {
            if (affine_.portFaceIds[static_cast<std::size_t>(i)] == faceId) {
                ++modeIndex;
            }
        }
        return multi->modes[static_cast<std::size_t>(modeIndex)];
    }
    return portModeSolver_.solve(faceId);
}

}  // namespace fem::sweep
