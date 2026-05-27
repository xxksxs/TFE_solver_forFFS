#include "bpfem/linalg/BiCGStabSolver.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace fem {

void BiCGStabSolver::setPreconditioner(std::shared_ptr<linalg::IPreconditioner> precon) {
    precon_ = std::move(precon);
}

SolveResult BiCGStabSolver::solve(const SparseMatrix& a, const std::vector<std::complex<double>>& b, int maxIterations, double tolerance) const {
    // Right-preconditioned BiCGSTAB (Saad 2003, Alg. 7.7) with optional
    // identity preconditioner. The classical recurrence is
    //   p_k = r_{k-1} + beta * (p_{k-1} - omega_{k-1} v_{k-1})
    //   v_k = A * M^{-1} p_k
    //   alpha_k = rho_k / <r0_tilde, v_k>
    //   s_k = r_{k-1} - alpha_k v_k
    //   t_k = A * M^{-1} s_k
    //   omega_k = <t_k, s_k> / <t_k, t_k>
    //   x_k = x_{k-1} + alpha_k M^{-1} p_k + omega_k M^{-1} s_k
    //   r_k = s_k - omega_k t_k
    // Below `pUse` is M^{-1} p_k and `sUse` is M^{-1} s_k. r0 (the shadow
    // residual) is fixed at r_0 = b - A x_0 = b for x_0 = 0.
    const std::size_t n = a.size();
    std::vector<std::complex<double>> x(n, 0.0), r = b, r0 = r, p(n, 0.0), v(n, 0.0), s(n), t(n);
    // p_hat = M^{-1} p, s_hat = M^{-1} s. When no preconditioner is set, we
    // skip the apply() call entirely and these aliases share storage with
    // p / s, so the algorithm collapses to the original unpreconditioned
    // form bit-for-bit.
    std::vector<std::complex<double>> pHat, sHat;
    if (precon_) {
        precon_->setup(a);
        pHat.assign(n, 0.0);
        sHat.assign(n, 0.0);
    }
    auto applyMinv = [&](const std::vector<std::complex<double>>& in,
                         std::vector<std::complex<double>>& out) -> const std::vector<std::complex<double>>& {
        if (!precon_) {
            return in;  // identity preconditioner; no copy
        }
        precon_->apply(in, out);
        return out;
    };

    std::complex<double> rhoPrev = 1.0, alpha = 1.0, omega = 1.0;
    const double normB = std::max(norm(b), 1.0e-30);
    double residual = norm(r) / normB;

    int iter = 0;
    for (; iter < maxIterations && residual > tolerance; ++iter) {
        // Standard BiCGSTAB recurrence; comments above mirror this loop.
        const std::complex<double> rho = dot(r0, r);
        if (std::abs(rho) < 1.0e-30) {
            break;  // breakdown: shadow residual lost the original direction
        }
        const std::complex<double> beta = (rho / rhoPrev) * (alpha / omega);
        for (std::size_t i = 0; i < n; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        const auto& pUse = applyMinv(p, pHat);
        v = a.multiply(pUse);
        const std::complex<double> denom = dot(r0, v);
        if (std::abs(denom) < 1.0e-30) {
            break;  // breakdown
        }
        alpha = rho / denom;
        for (std::size_t i = 0; i < n; ++i) {
            s[i] = r[i] - alpha * v[i];
        }
        // Early termination: BiCGSTAB's classical "lucky breakout" --
        // ||s|| already small, just take a half step and stop.
        if (norm(s) / normB < tolerance) {
            for (std::size_t i = 0; i < n; ++i) {
                x[i] += alpha * pUse[i];
            }
            r = s;
            break;
        }
        const auto& sUse = applyMinv(s, sHat);
        t = a.multiply(sUse);
        const std::complex<double> tt = dot(t, t);
        if (std::abs(tt) < 1.0e-30) {
            break;
        }
        omega = dot(t, s) / tt;
        for (std::size_t i = 0; i < n; ++i) {
            x[i] += alpha * pUse[i] + omega * sUse[i];
            r[i] = s[i] - omega * t[i];
        }
        residual = norm(r) / normB;
        rhoPrev = rho;
        if (std::abs(omega) < 1.0e-30) {
            break;  // stagnation
        }
    }
    return SolveResult{x, iter, residual};
}

std::complex<double> BiCGStabSolver::dot(const std::vector<std::complex<double>>& a, const std::vector<std::complex<double>>& b) {
    std::complex<double> value = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        value += std::conj(a[i]) * b[i];
    }
    return value;
}

double BiCGStabSolver::norm(const std::vector<std::complex<double>>& x) {
    double value = 0.0;
    for (const auto& xi : x) {
        value += std::norm(xi);
    }
    return std::sqrt(value);
}

}  // namespace fem

