#include "bpfem/linalg/BiCGStabSolver.hpp"
#include "bpfem/linalg/SolveResidual.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace fem {

namespace {

bool isFinite(std::complex<double> value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

}  // namespace

void BiCGStabSolver::setPreconditioner(std::shared_ptr<linalg::IPreconditioner> precon) {
    if (preparedMatrix_ != nullptr) {
        throw std::logic_error("Cannot replace the preconditioner in a prepared solve");
    }
    precon_ = std::move(precon);
}

void BiCGStabSolver::beginMatrix(const SparseMatrix& matrix) {
    if (preparedMatrix_ != nullptr) {
        throw std::logic_error("Iterative solver already has a prepared matrix");
    }
    auto matrixOperator = std::make_unique<linalg::IterativeOperator>(matrix, useMklSpmv_);
    if (precon_) {
        precon_->setup(matrix);
    }
    preparedMatrix_ = &matrix;
    preparedOperator_ = std::move(matrixOperator);
}

void BiCGStabSolver::endMatrix() noexcept {
    if (precon_) {
        precon_->clear();
    }
    preparedOperator_.reset();
    preparedMatrix_ = nullptr;
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
    // residual) is fixed at b / ||b|| for x_0 = 0.
    const std::size_t n = a.size();
    linalg::validateIterativeInput(a, b, maxIterations, tolerance);
    if (preparedMatrix_ != nullptr && preparedMatrix_ != &a) {
        throw std::logic_error("BiCGSTAB solve does not match its prepared matrix");
    }
    const double normB = linalg::vectorNorm(b);
    SolveResult result;
    auto& x = result.field;
    x.resize(n);
    if (normB == 0.0 || !std::isfinite(normB)) {
        linalg::checkResidual(a, b, tolerance, result);
        return result;
    }
    const linalg::IterativeOperator singleSolveOperator(a, useMklSpmv_ && !preparedOperator_);
    const auto& matrixOperator = preparedOperator_ ? *preparedOperator_ : singleSolveOperator;
    // Solve for x / ||b|| so inner products remain usable when moment RHS
    // magnitudes vary by many orders. Check the rescaled field on A x = b.
    std::vector<std::complex<double>> r = b;
    for (auto& value : r) {
        value /= normB;
    }
    const auto r0 = r;
    std::vector<std::complex<double>> p(n), v(n), s(n), t(n);
    // p_hat = M^{-1} p, s_hat = M^{-1} s. When no preconditioner is set, we
    // skip the apply() call entirely and these aliases share storage with
    // p / s.
    std::vector<std::complex<double>> pHat, sHat;
    if (precon_) {
        if (preparedMatrix_ == nullptr) {
            precon_->setup(a);
        }
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
    double residual = 1.0;

    int& iter = result.iterations;
    SolveStatus failure = SolveStatus::IterationLimit;
    while (iter < maxIterations && std::isfinite(residual) && residual > tolerance) {
        ++iter;
        // Standard BiCGSTAB recurrence; comments above mirror this loop.
        const std::complex<double> rho = dot(r0, r);
        if (!isFinite(rho) || rho == 0.0) {
            failure = isFinite(rho) ? SolveStatus::Breakdown : SolveStatus::NonFinite;
            break;  // breakdown: shadow residual lost the original direction
        }
        const std::complex<double> beta = (rho / rhoPrev) * (alpha / omega);
        if (!isFinite(beta)) {
            failure = SolveStatus::NonFinite;
            break;
        }
        for (std::size_t i = 0; i < n; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        const auto& pUse = applyMinv(p, pHat);
        v = matrixOperator.multiply(pUse);
        const std::complex<double> denom = dot(r0, v);
        if (!isFinite(denom) || denom == 0.0) {
            failure = isFinite(denom) ? SolveStatus::Breakdown : SolveStatus::NonFinite;
            break;  // breakdown
        }
        alpha = rho / denom;
        if (!isFinite(alpha)) {
            failure = SolveStatus::NonFinite;
            break;
        }
        for (std::size_t i = 0; i < n; ++i) {
            s[i] = r[i] - alpha * v[i];
        }
        // Early termination: BiCGSTAB's classical "lucky breakout" --
        // ||s|| already small, just take a half step and stop.
        if (linalg::vectorNorm(s) <= tolerance) {
            for (std::size_t i = 0; i < n; ++i) {
                x[i] += alpha * pUse[i];
            }
            failure = SolveStatus::ResidualTooLarge;
            break;
        }
        const auto& sUse = applyMinv(s, sHat);
        t = matrixOperator.multiply(sUse);
        const std::complex<double> tt = dot(t, t);
        if (!isFinite(tt) || tt == 0.0) {
            failure = isFinite(tt) ? SolveStatus::Breakdown : SolveStatus::NonFinite;
            break;
        }
        omega = dot(t, s) / tt;
        if (!isFinite(omega)) {
            failure = SolveStatus::NonFinite;
            break;
        }
        for (std::size_t i = 0; i < n; ++i) {
            x[i] += alpha * pUse[i] + omega * sUse[i];
            r[i] = s[i] - omega * t[i];
        }
        residual = linalg::vectorNorm(r);
        if (!std::isfinite(residual)) {
            failure = SolveStatus::NonFinite;
            break;
        }
        if (residual <= tolerance) {
            failure = SolveStatus::ResidualTooLarge;
            break;
        }
        rhoPrev = rho;
        if (omega == 0.0) {
            failure = SolveStatus::Breakdown;
            break;  // stagnation
        }
    }
    for (auto& value : x) {
        value *= normB;
    }
    linalg::checkResidual(a, b, tolerance, result, failure);
    return result;
}

std::complex<double> BiCGStabSolver::dot(const std::vector<std::complex<double>>& a, const std::vector<std::complex<double>>& b) {
    std::complex<double> value = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        value += std::conj(a[i]) * b[i];
    }
    return value;
}

}  // namespace fem

