#include "bpfem/linalg/IterativeOperator.hpp"

#include <stdexcept>

#ifdef BPFEM_USE_MKL
#include <mkl.h>
#endif

namespace fem::linalg {

struct IterativeOperator::State {
#ifdef BPFEM_USE_MKL
    sparse_matrix_t handle = nullptr;
    std::vector<MKL_INT> offsets;
    std::vector<MKL_INT> columns;
    std::vector<MKL_Complex16> values;
    std::vector<MKL_Complex16> input;
    std::vector<MKL_Complex16> output;
    matrix_descr descriptor{SPARSE_MATRIX_TYPE_SYMMETRIC, SPARSE_FILL_MODE_UPPER,
                            SPARSE_DIAG_NON_UNIT};

    ~State() {
        if (handle != nullptr) {
            mkl_sparse_destroy(handle);
        }
    }
#endif
};

IterativeOperator::IterativeOperator(const SparseMatrix& matrix, bool useMkl)
    : matrix_(matrix) {
    if (!useMkl || matrix.size() == 0 || matrix.values().empty()) {
        return;
    }
#ifdef BPFEM_USE_MKL
    auto state = std::make_unique<State>();
    const auto& offsets = matrix.rowOffsets();
    const auto& columns = matrix.columnIndices();
    state->offsets.assign(offsets.begin(), offsets.end());
    state->columns.assign(columns.begin(), columns.end());
    state->values.reserve(matrix.values().size());
    for (const auto& value : matrix.values()) {
        state->values.push_back({value.real(), value.imag()});
    }
    state->input.resize(matrix.size());
    state->output.resize(matrix.size());
    const auto n = static_cast<MKL_INT>(matrix.size());
    const auto status = mkl_sparse_z_create_csr(
        &state->handle, SPARSE_INDEX_BASE_ZERO, n, n, state->offsets.data(),
        state->offsets.data() + 1, state->columns.data(), state->values.data());
    if (status != SPARSE_STATUS_SUCCESS) {
        throw std::runtime_error("oneMKL could not create a complex CSR operator");
    }
    state_ = std::move(state);
#else
    throw std::runtime_error("Iterative MKL matrix multiplication requires Intel oneMKL");
#endif
}

IterativeOperator::~IterativeOperator() = default;

std::vector<std::complex<double>> IterativeOperator::multiply(
    const std::vector<std::complex<double>>& input) const {
#ifdef BPFEM_USE_MKL
    if (state_) {
        if (input.size() != matrix_.size()) {
            throw std::invalid_argument("Iterative operator vector size mismatch");
        }
        for (std::size_t i = 0; i < input.size(); ++i) {
            state_->input[i] = {input[i].real(), input[i].imag()};
        }
        const auto status = mkl_sparse_z_mv(
            SPARSE_OPERATION_NON_TRANSPOSE, {1.0, 0.0}, state_->handle, state_->descriptor,
            state_->input.data(), {0.0, 0.0}, state_->output.data());
        if (status != SPARSE_STATUS_SUCCESS) {
            throw std::runtime_error("oneMKL complex sparse matrix multiplication failed");
        }
        std::vector<std::complex<double>> result;
        result.reserve(input.size());
        for (const auto& value : state_->output) {
            result.emplace_back(value.real, value.imag);
        }
        return result;
    }
#endif
    return matrix_.multiply(input);
}

}  // namespace fem::linalg
