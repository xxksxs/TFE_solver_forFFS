#include "bpfem/linalg/SparseMatrix.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace fem {

SparseMatrix::SparseMatrix(std::size_t n) : n_(n), assemblyRows_(n) {}

SparseMatrix::SparseMatrix(std::shared_ptr<const SparsePattern> pattern)
    : n_(pattern ? pattern->n : 0), pattern_(std::move(pattern)), compressed_(true) {
    if (!pattern_) {
        throw std::runtime_error("SparseMatrix fixed pattern is null");
    }
    values_.assign(pattern_->columnIndices.size(), 0.0);
}

std::shared_ptr<const SparsePattern> SparseMatrix::createPattern(std::size_t n, const std::vector<std::vector<int>>& columnsByRow) {
    if (columnsByRow.size() != n) {
        throw std::runtime_error("SparseMatrix::createPattern row count mismatch");
    }
    auto pattern = std::make_shared<SparsePattern>();
    pattern->n = n;
    pattern->rowOffsets.assign(n + 1, 0);
    pattern->entryIndexByRow.resize(n);
    for (std::size_t r = 0; r < n; ++r) {
        pattern->rowOffsets[r] = static_cast<int>(pattern->columnIndices.size());
        int previous = -1;
        for (int col : columnsByRow[r]) {
            if (col < static_cast<int>(r) || static_cast<std::size_t>(col) >= n) {
                throw std::runtime_error("SparseMatrix::createPattern column out of upper-triangular range");
            }
            if (col == previous) {
                continue;
            }
            const std::size_t index = pattern->columnIndices.size();
            pattern->columnIndices.push_back(col);
            pattern->entryIndexByRow[r][col] = index;
            previous = col;
        }
    }
    pattern->rowOffsets[n] = static_cast<int>(pattern->columnIndices.size());
    return pattern;
}

std::size_t SparseMatrix::size() const {
    return n_;
}

void SparseMatrix::add(int r, int c, const std::complex<double>& v) {
    if (r < 0 || c < 0 || static_cast<std::size_t>(r) >= n_ || static_cast<std::size_t>(c) >= n_) {
        throw std::out_of_range("SparseMatrix::add index out of range");
    }
    if (c < r) {
        std::swap(r, c);
    }
    if (pattern_) {
        const auto& rowIndex = pattern_->entryIndexByRow[static_cast<std::size_t>(r)];
        const auto it = rowIndex.find(c);
        if (it == rowIndex.end()) {
            throw std::runtime_error("SparseMatrix::add entry is outside the prebuilt sparsity pattern");
        }
        values_[it->second] += v;
        return;
    }
    assemblyRows_[static_cast<std::size_t>(r)][c] += v;
    compressed_ = false;
}

void SparseMatrix::imposeDirichlet(int row, const std::complex<double>& value, std::vector<std::complex<double>>& rhs) {
    if (pattern_) {
        if (value != std::complex<double>(0.0, 0.0)) {
            throw std::runtime_error("SparseMatrix::imposeDirichlet with a nonzero value is not supported for fixed patterns");
        }
        imposeZeroDirichlet(std::unordered_set<int>{row}, rhs);
        return;
    }
    if (row < 0 || static_cast<std::size_t>(row) >= n_) {
        throw std::out_of_range("SparseMatrix::imposeDirichlet index out of range");
    }
    for (std::size_t r = 0; r < n_; ++r) {
        if (static_cast<int>(r) == row) {
            continue;
        }
        int storedRow = static_cast<int>(r);
        int storedCol = row;
        if (storedCol < storedRow) {
            std::swap(storedRow, storedCol);
        }
        auto& values = assemblyRows_[static_cast<std::size_t>(storedRow)];
        const auto it = values.find(storedCol);
        if (it != values.end()) {
            rhs[r] -= it->second * value;
            values.erase(it);
        }
    }
    assemblyRows_[static_cast<std::size_t>(row)].clear();
    for (std::size_t r = 0; r < static_cast<std::size_t>(row); ++r) {
        assemblyRows_[r].erase(row);
    }
    assemblyRows_[static_cast<std::size_t>(row)][row] = 1.0;
    rhs[static_cast<std::size_t>(row)] = value;
    compressed_ = false;
}

void SparseMatrix::imposeZeroDirichlet(const std::unordered_set<int>& rows, std::vector<std::complex<double>>& rhs) {
    for (int row : rows) {
        if (row < 0 || static_cast<std::size_t>(row) >= n_) {
            throw std::out_of_range("SparseMatrix::imposeZeroDirichlet index out of range");
        }
    }
    if (pattern_) {
        for (int row : rows) {
            const auto rowIndex = static_cast<std::size_t>(row);
            const int begin = pattern_->rowOffsets[rowIndex];
            const int end = pattern_->rowOffsets[rowIndex + 1];
            for (int k = begin; k < end; ++k) {
                values_[static_cast<std::size_t>(k)] = 0.0;
            }
            const auto it = pattern_->entryIndexByRow[rowIndex].find(row);
            if (it == pattern_->entryIndexByRow[rowIndex].end()) {
                throw std::runtime_error("SparseMatrix::imposeZeroDirichlet missing diagonal in fixed pattern");
            }
            values_[it->second] = 1.0;
            rhs[rowIndex] = 0.0;
        }
        return;
    }
    for (std::size_t r = 0; r < n_; ++r) {
        auto& values = assemblyRows_[r];
        if (rows.count(static_cast<int>(r)) != 0) {
            values.clear();
            continue;
        }
        for (auto it = values.begin(); it != values.end();) {
            if (rows.count(it->first) != 0) {
                it = values.erase(it);
            } else {
                ++it;
            }
        }
    }
    for (int row : rows) {
        assemblyRows_[static_cast<std::size_t>(row)][row] = 1.0;
        rhs[static_cast<std::size_t>(row)] = 0.0;
    }
    compressed_ = false;
}

std::vector<std::complex<double>> SparseMatrix::multiply(const std::vector<std::complex<double>>& x) const {
    if (x.size() != n_) {
        throw std::runtime_error("SparseMatrix::multiply vector size mismatch");
    }
    compress();
    const auto& offsets = rowOffsets();
    const auto& columns = columnIndices();
    std::vector<std::complex<double>> y(n_, 0.0);
    for (std::size_t r = 0; r < n_; ++r) {
        const int begin = offsets[r];
        const int end = offsets[r + 1];
        for (int k = begin; k < end; ++k) {
            const int c = columns[static_cast<std::size_t>(k)];
            const auto value = values_[static_cast<std::size_t>(k)];
            y[r] += value * x[static_cast<std::size_t>(c)];
            if (c != static_cast<int>(r)) {
                y[static_cast<std::size_t>(c)] += value * x[r];
            }
        }
    }
    return y;
}

void SparseMatrix::compress() const {
    if (pattern_) {
        return;
    }
    if (compressed_) {
        return;
    }
    rowOffsets_.assign(n_ + 1, 0);
    columnIndices_.clear();
    values_.clear();
    for (std::size_t r = 0; r < n_; ++r) {
        rowOffsets_[r] = static_cast<int>(values_.size());
        std::vector<std::pair<int, std::complex<double>>> entries(assemblyRows_[r].begin(), assemblyRows_[r].end());
        std::sort(entries.begin(), entries.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.first < rhs.first;
        });
        for (const auto& [col, value] : entries) {
            if (std::abs(value) <= 0.0) {
                continue;
            }
            columnIndices_.push_back(col);
            values_.push_back(value);
        }
    }
    rowOffsets_[n_] = static_cast<int>(values_.size());
    compressed_ = true;
}

bool SparseMatrix::isCompressed() const {
    return compressed_;
}

const std::vector<int>& SparseMatrix::rowOffsets() const {
    if (pattern_) {
        return pattern_->rowOffsets;
    }
    compress();
    return rowOffsets_;
}

const std::vector<int>& SparseMatrix::columnIndices() const {
    if (pattern_) {
        return pattern_->columnIndices;
    }
    compress();
    return columnIndices_;
}

const std::vector<std::complex<double>>& SparseMatrix::values() const {
    compress();
    return values_;
}

}  // namespace fem
