#include "bpfem/linalg/PreconSchwarz.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace fem::linalg {
namespace {

// Cuthill-McKee ordering keeps each contiguous core close in the matrix graph.
std::vector<int> orderGraph(const std::vector<std::vector<int>>& neighbors) {
    std::vector<int> roots(neighbors.size());
    std::iota(roots.begin(), roots.end(), 0);
    const auto lessDegree = [&](int a, int b) {
        const auto da = neighbors[a].size();
        const auto db = neighbors[b].size();
        return da == db ? a < b : da < db;
    };
    std::sort(roots.begin(), roots.end(), lessDegree);
    std::vector<int> order;
    order.reserve(roots.size());
    std::vector<bool> seen(roots.size());
    std::vector<int> next;
    for (const int root : roots) {
        if (seen[root]) {
            continue;
        }
        std::size_t cursor = order.size();
        order.push_back(root);
        seen[root] = true;
        while (cursor < order.size()) {
            next.clear();
            for (const int neighbor : neighbors[order[cursor++]]) {
                if (!seen[neighbor]) {
                    seen[neighbor] = true;
                    next.push_back(neighbor);
                }
            }
            std::sort(next.begin(), next.end(), lessDegree);
            order.insert(order.end(), next.begin(), next.end());
        }
    }
    return order;
}

}  // namespace

PreconSchwarz::PreconSchwarz(int subdomains, int overlap, LocalSolverFactory createLocalSolver)
    : subdomainCount_(subdomains), overlap_(overlap),
      createLocalSolver_(std::move(createLocalSolver)) {
    if (subdomains < 1 || overlap < 0 || !createLocalSolver_) {
        throw std::invalid_argument("Invalid Schwarz preconditioner configuration");
    }
}

void PreconSchwarz::clear() noexcept {
    subdomains_.clear();
    dimension_ = 0;
}

void PreconSchwarz::setup(const SparseMatrix& matrix) {
    clear();
    const auto n = matrix.size();
    if (n == 0) {
        return;
    }
    const auto& offsets = matrix.rowOffsets();
    const auto& columns = matrix.columnIndices();
    const auto& values = matrix.values();
    std::vector<std::vector<int>> neighbors(n);
    for (std::size_t row = 0; row < n; ++row) {
        for (int k = offsets[row]; k < offsets[row + 1]; ++k) {
            const auto col = columns[k];
            if (static_cast<std::size_t>(col) != row && values[k] != 0.0) {
                neighbors[row].push_back(col);
                neighbors[col].push_back(static_cast<int>(row));
            }
        }
    }
    const auto order = orderGraph(neighbors);
    const auto count = std::min(n, static_cast<std::size_t>(subdomainCount_));
    std::vector<Subdomain> domains(count);
    std::vector<int> localIndex(n, -1);
    for (std::size_t part = 0; part < count; ++part) {
        auto& domain = domains[part];
        // Split evenly without computing n * part.
        const auto start = (n / count) * part + std::min(part, n % count);
        domain.coreSize = n / count + (part < n % count ? 1 : 0);
        domain.indices.assign(order.begin() + start, order.begin() + start + domain.coreSize);
        for (std::size_t i = 0; i < domain.coreSize; ++i) {
            localIndex[domain.indices[i]] = static_cast<int>(i);
        }
        std::size_t layerStart = 0;
        for (int layer = 0; layer < overlap_ && layerStart < domain.indices.size(); ++layer) {
            const auto layerEnd = domain.indices.size();
            for (std::size_t i = layerStart; i < layerEnd; ++i) {
                for (const int neighbor : neighbors[domain.indices[i]]) {
                    if (localIndex[neighbor] < 0) {
                        localIndex[neighbor] = static_cast<int>(domain.indices.size());
                        domain.indices.push_back(neighbor);
                    }
                }
            }
            layerStart = layerEnd;
        }
        SparseMatrix local(domain.indices.size());
        for (std::size_t i = 0; i < domain.indices.size(); ++i) {
            const auto global = domain.indices[i];
            for (int k = offsets[global]; k < offsets[global + 1]; ++k) {
                const auto j = localIndex[columns[k]];
                if (j >= 0) {
                    local.add(static_cast<int>(i), j, values[k]);
                }
            }
        }
        domain.solver = createLocalSolver_();
        try {
            domain.solver->factorize(local);
        } catch (const std::exception& error) {
            throw std::runtime_error("Schwarz subdomain " + std::to_string(part + 1)
                + " factorization failed: " + error.what());
        }
        for (const int global : domain.indices) {
            localIndex[global] = -1;
        }
    }
    dimension_ = n;
    subdomains_ = std::move(domains);
}

void PreconSchwarz::apply(const std::vector<std::complex<double>>& x,
                         std::vector<std::complex<double>>& y) const {
    if (x.size() != dimension_) {
        throw std::invalid_argument("Schwarz input size does not match the prepared matrix");
    }
    std::vector<std::complex<double>> original;
    if (&x == &y) {
        original = x;
    }
    const auto& source = &x == &y ? original : x;
    y.resize(dimension_);
    std::vector<std::complex<double>> rhs;
    for (const auto& domain : subdomains_) {
        rhs.resize(domain.indices.size());
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            rhs[i] = source[domain.indices[i]];
        }
        const auto solution = domain.solver->solveFactorized(rhs);
        for (std::size_t i = 0; i < domain.coreSize; ++i) {
            y[domain.indices[i]] = solution.field[i];
        }
    }
}

}  // namespace fem::linalg
