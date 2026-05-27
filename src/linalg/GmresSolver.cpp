#include "bpfem/linalg/GmresSolver.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace fem {

namespace {

using Complex = std::complex<double>;

// Hermitian inner product for orthonormalization.
Complex hdot(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    Complex s(0.0, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        s += std::conj(a[i]) * b[i];
    }
    return s;
}

double norm2(const std::vector<Complex>& v) {
    double s = 0.0;
    for (const auto& z : v) {
        s += std::norm(z);
    }
    return std::sqrt(s);
}

// Compute a Givens rotation that zeros y when applied to (x, y)^T:
// [c   s; -conj(s)  c] (x; y) = (r; 0).
// Standard Bindel/Demmel BLAS ZROTG-style formula stable for complex inputs.
void givens(Complex x, Complex y, double& c, Complex& s, Complex& r) {
    const double ax = std::abs(x);
    const double ay = std::abs(y);
    if (ay == 0.0) {
        c = 1.0;
        s = Complex(0.0, 0.0);
        r = x;
        return;
    }
    if (ax == 0.0) {
        c = 0.0;
        // s such that |s|=1; choose s = conj(y)/|y|
        s = std::conj(y) / ay;
        r = Complex(ay, 0.0);
        return;
    }
    const double rho = std::sqrt(ax * ax + ay * ay);
    c = ax / rho;
    // s = (x / |x|) * conj(y) / rho
    s = (x / ax) * std::conj(y) / rho;
    r = (x / ax) * rho;
}

void applyGivens(double c, Complex s, Complex& a, Complex& b) {
    const Complex newA = c * a + s * b;
    const Complex newB = -std::conj(s) * a + c * b;
    a = newA;
    b = newB;
}

}  // namespace

void GmresSolver::setPreconditioner(std::shared_ptr<linalg::IPreconditioner> precon) {
    precon_ = std::move(precon);
}

SolveResult GmresSolver::solve(const SparseMatrix& A,
                               const std::vector<Complex>& b,
                               int maxIterations,
                               double tolerance) const {
    const std::size_t n = A.size();
    const int m = restartLength_;

    std::vector<Complex> x(n, Complex(0.0, 0.0));
    std::vector<Complex> r(n);

    if (precon_) {
        precon_->setup(A);
    }

    auto applyPrecon = [&](const std::vector<Complex>& in,
                            std::vector<Complex>& scratch) -> const std::vector<Complex>& {
        if (!precon_) {
            return in;
        }
        precon_->apply(in, scratch);
        return scratch;
    };

    const double normB = std::max(norm2(b), 1.0e-30);

    // Initial residual r0 = b - A x0 = b (x0 = 0).
    r = b;
    double residual = norm2(r) / normB;

    // Krylov basis V (m+1 columns, each length n) and Hessenberg matrix H
    // (column-major (m+1) x m). Givens rotations stored as (c, s) per column.
    std::vector<std::vector<Complex>> V;
    V.reserve(static_cast<std::size_t>(m + 1));
    std::vector<Complex> H(static_cast<std::size_t>((m + 1) * m), Complex(0.0, 0.0));
    std::vector<double> givensC(static_cast<std::size_t>(m), 0.0);
    std::vector<Complex> givensS(static_cast<std::size_t>(m), Complex(0.0, 0.0));
    std::vector<Complex> g(static_cast<std::size_t>(m + 1), Complex(0.0, 0.0));
    std::vector<Complex> preconScratch(n, Complex(0.0, 0.0));

    int iter = 0;
    while (iter < maxIterations && residual > tolerance) {
        // Restart: compute current residual.
        const auto Ax = A.multiply(x);
        for (std::size_t i = 0; i < n; ++i) {
            r[i] = b[i] - Ax[i];
        }
        const double beta = norm2(r);
        if (beta == 0.0) break;
        residual = beta / normB;
        if (residual <= tolerance) break;

        V.clear();
        V.emplace_back(n, Complex(0.0, 0.0));
        for (std::size_t i = 0; i < n; ++i) {
            V[0][i] = r[i] / beta;
        }
        std::fill(H.begin(), H.end(), Complex(0.0, 0.0));
        std::fill(givensC.begin(), givensC.end(), 0.0);
        std::fill(givensS.begin(), givensS.end(), Complex(0.0, 0.0));
        std::fill(g.begin(), g.end(), Complex(0.0, 0.0));
        g[0] = Complex(beta, 0.0);

        int kFinal = -1;
        for (int j = 0; j < m && iter < maxIterations; ++j) {
            ++iter;

            // w = A * M^{-1} v_j  (right preconditioning).
            const auto& vjPrec = applyPrecon(V[static_cast<std::size_t>(j)], preconScratch);
            std::vector<Complex> w = A.multiply(vjPrec);

            // Modified Gram-Schmidt with one pass of re-orthogonalization.
            // Both passes accumulate into H(i, j); the second pass refines
            // already-orthogonalized w against round-off (DGKS-style).
            for (int pass = 0; pass < 2; ++pass) {
                for (int i = 0; i <= j; ++i) {
                    const Complex h = hdot(V[static_cast<std::size_t>(i)], w);
                    H[static_cast<std::size_t>(i) * static_cast<std::size_t>(m)
                      + static_cast<std::size_t>(j)] += h;
                    for (std::size_t k = 0; k < n; ++k) {
                        w[k] -= h * V[static_cast<std::size_t>(i)][k];
                    }
                }
            }
            const double hnext = norm2(w);
            H[static_cast<std::size_t>(j + 1) * static_cast<std::size_t>(m)
              + static_cast<std::size_t>(j)] = Complex(hnext, 0.0);

            // Apply previous Givens rotations to column j.
            for (int i = 0; i < j; ++i) {
                Complex& a = H[static_cast<std::size_t>(i) * static_cast<std::size_t>(m)
                                + static_cast<std::size_t>(j)];
                Complex& bb = H[static_cast<std::size_t>(i + 1) * static_cast<std::size_t>(m)
                                + static_cast<std::size_t>(j)];
                applyGivens(givensC[static_cast<std::size_t>(i)],
                            givensS[static_cast<std::size_t>(i)], a, bb);
            }

            // New Givens rotation to zero H(j+1, j).
            Complex hjj = H[static_cast<std::size_t>(j) * static_cast<std::size_t>(m)
                             + static_cast<std::size_t>(j)];
            Complex hj1j = H[static_cast<std::size_t>(j + 1) * static_cast<std::size_t>(m)
                              + static_cast<std::size_t>(j)];
            double c = 1.0;
            Complex s(0.0, 0.0), rPiv(0.0, 0.0);
            givens(hjj, hj1j, c, s, rPiv);
            givensC[static_cast<std::size_t>(j)] = c;
            givensS[static_cast<std::size_t>(j)] = s;
            H[static_cast<std::size_t>(j) * static_cast<std::size_t>(m)
              + static_cast<std::size_t>(j)] = rPiv;
            H[static_cast<std::size_t>(j + 1) * static_cast<std::size_t>(m)
              + static_cast<std::size_t>(j)] = Complex(0.0, 0.0);

            // Update RHS g.
            applyGivens(c, s, g[static_cast<std::size_t>(j)],
                        g[static_cast<std::size_t>(j + 1)]);

            residual = std::abs(g[static_cast<std::size_t>(j + 1)]) / normB;
            kFinal = j;

            if (residual <= tolerance) {
                break;
            }

            // Lucky breakdown / next basis vector.
            if (hnext == 0.0) {
                break;
            }
            std::vector<Complex> next(n);
            for (std::size_t k = 0; k < n; ++k) {
                next[k] = w[k] / hnext;
            }
            V.push_back(std::move(next));
        }

        if (kFinal < 0) {
            break;  // converged before any iteration produced an update
        }

        // Solve upper triangular H(0..kFinal, 0..kFinal) y = g(0..kFinal).
        const int kf = kFinal;
        std::vector<Complex> y(static_cast<std::size_t>(kf + 1), Complex(0.0, 0.0));
        for (int i = kf; i >= 0; --i) {
            Complex sum = g[static_cast<std::size_t>(i)];
            for (int k = i + 1; k <= kf; ++k) {
                sum -= H[static_cast<std::size_t>(i) * static_cast<std::size_t>(m)
                         + static_cast<std::size_t>(k)] * y[static_cast<std::size_t>(k)];
            }
            const Complex diag = H[static_cast<std::size_t>(i) * static_cast<std::size_t>(m)
                                    + static_cast<std::size_t>(i)];
            y[static_cast<std::size_t>(i)] = sum / diag;
        }

        // x_update = M^{-1} (V * y)
        std::vector<Complex> Vy(n, Complex(0.0, 0.0));
        for (int i = 0; i <= kf; ++i) {
            for (std::size_t k = 0; k < n; ++k) {
                Vy[k] += y[static_cast<std::size_t>(i)] * V[static_cast<std::size_t>(i)][k];
            }
        }
        const auto& update = applyPrecon(Vy, preconScratch);
        for (std::size_t k = 0; k < n; ++k) {
            x[k] += update[k];
        }
    }

    return SolveResult{x, iter, residual};
}

}  // namespace fem

