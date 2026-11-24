#pragma once

#include "bpfem/linalg/ISparseSolver.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace fem::linalg {

inline void validateIterativeInput(
    const SparseMatrix& matrix, const std::vector<std::complex<double>>& rhs,
    int maxIterations, double tolerance) {
    if (rhs.size() != matrix.size() || maxIterations < 1
        || !std::isfinite(tolerance) || tolerance <= 0.0) {
        throw std::invalid_argument("Invalid iterative solve dimensions or configuration");
    }
    for (const auto& value : rhs) {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
            throw std::invalid_argument("Non-finite iterative solve RHS");
        }
    }
}

class VectorNorm {
public:
    void add(std::complex<double> value) {
        addComponent(std::abs(value.real()));
        addComponent(std::abs(value.imag()));
    }

    double value() const { return scale_ * std::sqrt(squares_); }

private:
    void addComponent(double magnitude) {
        if (!std::isfinite(magnitude)) {
            scale_ = magnitude;
            return;
        }
        if (magnitude == 0.0 || !std::isfinite(scale_)) {
            return;
        }
        if (magnitude > scale_) {
            const double ratio = scale_ / magnitude;
            squares_ = 1.0 + squares_ * ratio * ratio;
            scale_ = magnitude;
        } else {
            const double ratio = magnitude / scale_;
            squares_ += ratio * ratio;
        }
    }

    double scale_ = 0.0;
    double squares_ = 0.0;
};

inline double vectorNorm(const std::vector<std::complex<double>>& vector) {
    VectorNorm norm;
    for (const auto value : vector) {
        norm.add(value);
    }
    return norm.value();
}

inline void checkResidual(
    const SparseMatrix& matrix, const std::vector<std::complex<double>>& rhs,
    double tolerance, SolveResult& result,
    SolveStatus failure = SolveStatus::ResidualTooLarge) {
    if (rhs.size() != matrix.size() || result.field.size() != matrix.size()) {
        throw std::invalid_argument("Linear solve returned incompatible vector dimensions");
    }
    const auto product = matrix.multiply(result.field);
    VectorNorm errorNorm;
    VectorNorm rhsNorm;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        errorNorm.add(product[i] - rhs[i]);
        rhsNorm.add(rhs[i]);
    }
    const double error = errorNorm.value();
    const double normB = rhsNorm.value();
    result.absoluteResidual = error;
    result.residual = normB == 0.0 ? error : error / normB;
    if (!std::isfinite(error) || !std::isfinite(normB) || !std::isfinite(result.residual)) {
        result.status = SolveStatus::NonFinite;
    } else if (result.residual <= tolerance) {
        result.status = SolveStatus::Converged;
    } else {
        result.status = failure;
    }
}

inline const char* solveStatusName(SolveStatus status) {
    switch (status) {
        case SolveStatus::Converged: return "converged";
        case SolveStatus::IterationLimit: return "iteration_limit";
        case SolveStatus::Breakdown: return "breakdown";
        case SolveStatus::NonFinite: return "non_finite";
        case SolveStatus::ResidualTooLarge: return "residual_too_large";
        case SolveStatus::Unchecked: return "unchecked";
    }
    return "unknown";
}

inline void requireConverged(const SolveResult& result, const char* backend) {
    if (result.status == SolveStatus::Converged) {
        return;
    }
    std::ostringstream message;
    message << backend << " linear solve failed (" << solveStatusName(result.status) << ", iterations="
            << result.iterations << ", relative residual=" << result.residual
            << ", absolute residual=" << result.absoluteResidual << ')';
    throw std::runtime_error(message.str());
}

}  // namespace fem::linalg
