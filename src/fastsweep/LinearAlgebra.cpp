#include "bpfem/fastsweep/LinearAlgebra.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace fem::fastsweep {

namespace {

// 将二维行列索引映射到一维行主序数组位置。
std::size_t idx(int row, int col, int n) {
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(n)
        + static_cast<std::size_t>(col);
}

}  // namespace

// 把稀疏实数端口向量提升为全阶复数向量，便于后续复数矩阵运算。
std::vector<Complex> liftRealSparseVector(const std::vector<std::pair<int, double>>& sparse,
                                          std::size_t n) {
    std::vector<Complex> out(n, Complex(0.0, 0.0));
    for (const auto& [i, w] : sparse) {
        if (i >= 0 && static_cast<std::size_t>(i) < n) {
            out[static_cast<std::size_t>(i)] = w;
        }
    }
    return out;
}

// 计算不取共轭的双线性内积，匹配复对称 FEM/Galerkin 投影约定。
Complex bilinear(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    Complex s(0.0, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        s += a[i] * b[i];
    }
    return s;
}

// 计算标准 Hermitian 内积，用于正交化、范数和稳定性诊断。
Complex hdot(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    Complex s(0.0, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        s += std::conj(a[i]) * b[i];
    }
    return s;
}

// 返回复向量的二范数。
double norm2(const std::vector<Complex>& a) {
    double s = 0.0;
    for (const auto& z : a) {
        s += std::norm(z);
    }
    return std::sqrt(s);
}

// 对候选向量执行改进 Gram-Schmidt 正交化，并返回剩余向量范数。
double modifiedGramSchmidt(std::vector<Complex>& w,
                           const std::vector<std::vector<Complex>>& basis,
                           int reorthogonalizationPasses) {
    const int passes = std::max(1, reorthogonalizationPasses);
    for (int pass = 0; pass < passes; ++pass) {
        for (const auto& v : basis) {
            const Complex h = hdot(v, w);
            for (std::size_t i = 0; i < w.size(); ++i) {
                w[i] -= h * v[i];
            }
        }
    }
    return norm2(w);
}

// 计算基向量集合相对单位正交关系的最大偏差。
double orthogonalityError(const std::vector<std::vector<Complex>>& basis) {
    double err = 0.0;
    for (std::size_t i = 0; i < basis.size(); ++i) {
        for (std::size_t j = 0; j < basis.size(); ++j) {
            const Complex expected = (i == j) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
            err = std::max(err, std::abs(hdot(basis[i], basis[j]) - expected));
        }
    }
    return err;
}

// 用带部分主元的稠密 LU 解小型复线性系统，供 ROM 和 Padé 小系统使用。
bool denseSolve(std::vector<Complex>& a,
                std::vector<Complex>& rhs,
                int n,
                double* pivotRatio) {
    double minPivot = std::numeric_limits<double>::infinity();
    double maxPivot = 0.0;

    for (int k = 0; k < n; ++k) {
        int pivotRow = k;
        double pivotMag = std::abs(a[idx(k, k, n)]);
        for (int r = k + 1; r < n; ++r) {
            const double mag = std::abs(a[idx(r, k, n)]);
            if (mag > pivotMag) {
                pivotMag = mag;
                pivotRow = r;
            }
        }
        if (pivotMag < 1.0e-30) {
            return false;
        }
        minPivot = std::min(minPivot, pivotMag);
        maxPivot = std::max(maxPivot, pivotMag);
        if (pivotRow != k) {
            for (int c = 0; c < n; ++c) {
                std::swap(a[idx(k, c, n)], a[idx(pivotRow, c, n)]);
            }
            std::swap(rhs[static_cast<std::size_t>(k)],
                      rhs[static_cast<std::size_t>(pivotRow)]);
        }

        const Complex pivot = a[idx(k, k, n)];
        for (int r = k + 1; r < n; ++r) {
            const Complex factor = a[idx(r, k, n)] / pivot;
            a[idx(r, k, n)] = factor;
            for (int c = k + 1; c < n; ++c) {
                a[idx(r, c, n)] -= factor * a[idx(k, c, n)];
            }
            rhs[static_cast<std::size_t>(r)] -= factor * rhs[static_cast<std::size_t>(k)];
        }
    }

    for (int r = n - 1; r >= 0; --r) {
        Complex sum = rhs[static_cast<std::size_t>(r)];
        for (int c = r + 1; c < n; ++c) {
            sum -= a[idx(r, c, n)] * rhs[static_cast<std::size_t>(c)];
        }
        rhs[static_cast<std::size_t>(r)] = sum / a[idx(r, r, n)];
    }

    if (pivotRatio != nullptr) {
        *pivotRatio = maxPivot / std::max(minPivot, 1.0e-300);
    }
    return true;
}

}  // namespace fem::fastsweep
