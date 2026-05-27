#include "bpfem/sweep/AlpsSweep.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/linalg/MklPardisoSolver.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fem::sweep {

namespace {

using Complex = std::complex<double>;

// Sparse symmetric mat-vec for matrices stored in upper-triangle CSR.
std::vector<Complex> spmv(const SparseMatrix& A, const std::vector<Complex>& x) {
    return A.multiply(x);
}

// dot using transpose (no conjugate) - the bilinear form for complex-symmetric
// Galerkin: V^T K V, where V^T is plain transpose.
Complex bilinear(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    Complex s(0.0, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        s += a[i] * b[i];
    }
    return s;
}

// Hermitian dot for orthonormalization stability.
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

// Modified Gram-Schmidt with one pass of reorthogonalization.
double mgs(std::vector<Complex>& w, const std::vector<std::vector<Complex>>& V, int numCols) {
    for (int pass = 0; pass < 2; ++pass) {
        for (int j = 0; j < numCols; ++j) {
            const Complex h = hdot(V[static_cast<std::size_t>(j)], w);
            for (std::size_t i = 0; i < w.size(); ++i) {
                w[i] -= h * V[static_cast<std::size_t>(j)][i];
            }
        }
    }
    return norm2(w);
}

// Lift a sparse (index, weight) port coupling vector to full DOF size.
std::vector<Complex> liftPortVector(const std::vector<std::pair<int, double>>& sparse, std::size_t n) {
    std::vector<Complex> out(n, Complex(0.0, 0.0));
    for (const auto& [i, w] : sparse) {
        if (i >= 0 && static_cast<std::size_t>(i) < n) {
            out[static_cast<std::size_t>(i)] = w;
        }
    }
    return out;
}

// Hand-rolled partial-pivoting LU for a small q x q dense complex system.
// q is at most a few hundred so this is good enough and stays MKL-LAPACK free.
bool denseSolve(std::vector<Complex>& A, std::vector<Complex>& rhs, int n) {
    auto idx = [n](int r, int c) {
        return static_cast<std::size_t>(r) * static_cast<std::size_t>(n) + static_cast<std::size_t>(c);
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
            std::swap(rhs[static_cast<std::size_t>(k)], rhs[static_cast<std::size_t>(pivotRow)]);
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
    for (int k = n - 1; k >= 0; --k) {
        Complex sum = rhs[static_cast<std::size_t>(k)];
        for (int c = k + 1; c < n; ++c) {
            sum -= A[idx(k, c)] * rhs[static_cast<std::size_t>(c)];
        }
        rhs[static_cast<std::size_t>(k)] = sum / A[idx(k, k)];
    }
    return true;
}

}  // namespace

AlpsSweep::AlpsSweep(const ProjectDefinition& project,
                     const FEMAssembler& assembler,
                     const PortModeSolver& portModeSolver,
                     AlpsOptions options)
    : project_(project), assembler_(assembler), portModeSolver_(portModeSolver), options_(options) {}

SweepResult AlpsSweep::run(const std::vector<double>& frequencies, const SweepContext& ctx) {
    const double expansionHz = options_.expansionFrequencyHz > 0.0
        ? options_.expansionFrequencyHz
        : (frequencies.empty()
            ? 0.5 * (project_.sweep.startHz + project_.sweep.endHz)
            : 0.5 * (frequencies.front() + frequencies.back()));

    ctx.log.info("ALPS expansion frequency: " + std::to_string(expansionHz / 1.0e9) + " GHz");
    ctx.log.info("ALPS Krylov order per port: " + std::to_string(options_.krylovOrder));
    const int romDim = buildOffline(expansionHz);
    ctx.log.info("ALPS ROM dimension: " + std::to_string(romDim));

    SweepResult out;
    out.points.reserve(frequencies.size());
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        const double f = frequencies[i];
        if (i == 0 || (i + 1) % 50 == 0 || i + 1 == frequencies.size()) {
            ctx.log.info("ALPS evaluating point " + std::to_string(i + 1) + "/"
                         + std::to_string(frequencies.size())
                         + " @ " + std::to_string(f / 1.0e9) + " GHz");
        }
        out.points.push_back(evaluate(f));
        out.lastFrequencyHz = f;
    }
    // ALPS does not retain a full-space basis: lastEdgeDofs is left empty.
    // Application interprets this as "no field VTU" and the SweepContext
    // onFieldSolved callback is intentionally not invoked.
    (void)ctx.onFieldSolved;
    return out;
}

int AlpsSweep::buildOffline(double expansionFrequencyHz) {
    expansionFrequencyHz_ = expansionFrequencyHz;
    numProjectPorts_ = static_cast<int>(project_.ports.size());
    if (numProjectPorts_ == 0) {
        throw std::runtime_error("AlpsSweep: no ports defined");
    }

    affine_ = assembler_.buildAffineSystem();
    fullDim_ = static_cast<int>(affine_.K.size());
    numVirtualPorts_ = static_cast<int>(affine_.portCoupling.size());
    if (numVirtualPorts_ == 0) {
        throw std::runtime_error("AlpsSweep: AffineSystem produced zero port entries");
    }
    if (!affine_.lossless) {
        throw std::runtime_error("AlpsSweep: lossy materials (sigma != 0) are not supported in this MVP. Use --sweep direct.");
    }

    std::vector<std::complex<double>> rhsTemp;
    SparseMatrix A0 = assembler_.assemble(expansionFrequencyHz_, rhsTemp);

    MklPardisoSolver directSolver;
    auto solveAt = [&](const std::vector<Complex>& b) -> std::vector<Complex> {
        SolveResult r = directSolver.solve(A0, b);
        return r.field;
    };

    std::vector<std::vector<Complex>> V;
    std::vector<std::vector<Complex>> portM;
    portM.reserve(static_cast<std::size_t>(numVirtualPorts_));
    for (int v = 0; v < numVirtualPorts_; ++v) {
        portM.push_back(liftPortVector(affine_.portCoupling[static_cast<std::size_t>(v)],
                                       static_cast<std::size_t>(fullDim_)));
    }

    const int targetCols = options_.krylovOrder * numVirtualPorts_;
    V.reserve(static_cast<std::size_t>(targetCols));

    auto appendIfIndependent = [&](std::vector<Complex>&& w) {
        const double pre = norm2(w);
        if (pre < options_.dropTolerance) {
            return;
        }
        const double post = mgs(w, V, static_cast<int>(V.size()));
        if (post < options_.dropTolerance * std::max(pre, 1.0)) {
            return;
        }
        for (auto& z : w) {
            z /= post;
        }
        V.push_back(std::move(w));
    };

    for (int v = 0; v < numVirtualPorts_; ++v) {
        std::vector<Complex> w = solveAt(portM[static_cast<std::size_t>(v)]);
        appendIfIndependent(std::move(w));
    }

    std::size_t cursor = 0;
    while (static_cast<int>(V.size()) < targetCols && cursor < V.size()) {
        const auto Mv = spmv(affine_.M, V[cursor]);
        std::vector<Complex> w = solveAt(Mv);
        appendIfIndependent(std::move(w));
        ++cursor;
    }

    romDim_ = static_cast<int>(V.size());
    if (romDim_ == 0) {
        throw std::runtime_error("AlpsSweep: Krylov subspace collapsed to zero dimension");
    }

    Ktilde_.assign(static_cast<std::size_t>(romDim_) * static_cast<std::size_t>(romDim_),
                   Complex(0.0, 0.0));
    Mtilde_.assign(static_cast<std::size_t>(romDim_) * static_cast<std::size_t>(romDim_),
                   Complex(0.0, 0.0));
    for (int j = 0; j < romDim_; ++j) {
        const auto Kvj = spmv(affine_.K, V[static_cast<std::size_t>(j)]);
        const auto Mvj = spmv(affine_.M, V[static_cast<std::size_t>(j)]);
        for (int i = 0; i < romDim_; ++i) {
            Ktilde_[static_cast<std::size_t>(i) * static_cast<std::size_t>(romDim_)
                    + static_cast<std::size_t>(j)] =
                bilinear(V[static_cast<std::size_t>(i)], Kvj);
            Mtilde_[static_cast<std::size_t>(i) * static_cast<std::size_t>(romDim_)
                    + static_cast<std::size_t>(j)] =
                bilinear(V[static_cast<std::size_t>(i)], Mvj);
        }
    }

    portModeReduced_.assign(static_cast<std::size_t>(numVirtualPorts_), {});
    for (int v = 0; v < numVirtualPorts_; ++v) {
        portModeReduced_[static_cast<std::size_t>(v)].assign(static_cast<std::size_t>(romDim_),
                                                              Complex(0.0, 0.0));
        for (int i = 0; i < romDim_; ++i) {
            portModeReduced_[static_cast<std::size_t>(v)][static_cast<std::size_t>(i)] =
                bilinear(V[static_cast<std::size_t>(i)], portM[static_cast<std::size_t>(v)]);
        }
    }

    ready_ = true;
    return romDim_;
}

SParameterPoint AlpsSweep::evaluate(double frequencyHz) const {
    if (!ready_) {
        throw std::runtime_error("AlpsSweep::evaluate called before buildOffline");
    }

    SParameterPoint sp;
    sp.frequencyHz = frequencyHz;
    if (numProjectPorts_ < 2) {
        return sp;
    }

    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double k0Sq = k0 * k0;

    std::vector<double> beta(static_cast<std::size_t>(numVirtualPorts_), 0.0);
    std::vector<double> sNorm(static_cast<std::size_t>(numVirtualPorts_), 0.0);
    std::vector<Complex> incident(static_cast<std::size_t>(numVirtualPorts_), Complex(0.0, 0.0));

    std::vector<int> dominantVirtualByProject(static_cast<std::size_t>(numProjectPorts_), -1);

    for (int v = 0; v < numVirtualPorts_; ++v) {
        const double kc2 = affine_.portCutoffSquared[static_cast<std::size_t>(v)];
        const double betaSq = k0Sq - kc2;
        beta[static_cast<std::size_t>(v)] = (betaSq > 0.0) ? std::sqrt(betaSq) : 0.0;

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

        if (affine_.isExcitationMode[static_cast<std::size_t>(v)]) {
            const auto& port = project_.ports[static_cast<std::size_t>(
                affine_.projectPortIndex[static_cast<std::size_t>(v)])];
            const double s = sNorm[static_cast<std::size_t>(v)];
            if (s > 0.0) {
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
        if (bv <= 0.0) {
            continue;
        }
        if (incident[static_cast<std::size_t>(v)] == Complex(0.0, 0.0)) {
            continue;
        }
        const Complex factor = Complex(0.0, 2.0) * bv * incident[static_cast<std::size_t>(v)];
        const auto& vp = portModeReduced_[static_cast<std::size_t>(v)];
        for (int i = 0; i < q; ++i) {
            rhs[static_cast<std::size_t>(i)] += factor * vp[static_cast<std::size_t>(i)];
        }
    }

    std::vector<Complex> AtildeWork = Atilde;
    std::vector<Complex> xTilde = rhs;
    if (!denseSolve(AtildeWork, xTilde, q)) {
        throw std::runtime_error("AlpsSweep::evaluate: reduced system singular at "
                                 + std::to_string(frequencyHz) + " Hz");
    }

    auto portProjection = [&](int virtualIdx) {
        const auto& vp = portModeReduced_[static_cast<std::size_t>(virtualIdx)];
        Complex sum(0.0, 0.0);
        for (int i = 0; i < q; ++i) {
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

}  // namespace fem::sweep

