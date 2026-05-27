#pragma once

#include <complex>
#include <cstddef>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fem {

struct SparsePattern {
    std::size_t n = 0;
    std::vector<int> rowOffsets;
    std::vector<int> columnIndices;
    std::vector<std::unordered_map<int, std::size_t>> entryIndexByRow;
};

struct SparseMatrix {
    explicit SparseMatrix(std::size_t n = 0);
    explicit SparseMatrix(std::shared_ptr<const SparsePattern> pattern);

    static std::shared_ptr<const SparsePattern> createPattern(std::size_t n, const std::vector<std::vector<int>>& columnsByRow);

    std::size_t size() const;
    void add(int r, int c, const std::complex<double>& v);
    void imposeDirichlet(int row, const std::complex<double>& value, std::vector<std::complex<double>>& rhs);
    void imposeZeroDirichlet(const std::unordered_set<int>& rows, std::vector<std::complex<double>>& rhs);
    std::vector<std::complex<double>> multiply(const std::vector<std::complex<double>>& x) const;
    void compress() const;
    bool isCompressed() const;
    const std::vector<int>& rowOffsets() const;
    const std::vector<int>& columnIndices() const;
    const std::vector<std::complex<double>>& values() const;

private:
    std::size_t n_ = 0;
    std::shared_ptr<const SparsePattern> pattern_;
    std::vector<std::unordered_map<int, std::complex<double>>> assemblyRows_;
    mutable bool compressed_ = false;
    mutable std::vector<int> rowOffsets_;
    mutable std::vector<int> columnIndices_;
    mutable std::vector<std::complex<double>> values_;
};

}  // namespace fem
