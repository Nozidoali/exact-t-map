#include "solver/parity-matrix.hpp"

#include <algorithm>
#include <functional>

namespace exact_t {

ParityMatrix::ParityMatrix(uint32_t n, uint32_t m) : n_(n), m_(m) {
    matrix_.resize(n, Bits(m + n));
    for (uint32_t i = 0; i < n; ++i)
        matrix_[i].set(m + i);
}

ParityMatrix::ParityMatrix(std::vector<std::pair<uint32_t, int>> const& terms, uint32_t n)
    : n_(n), m_(static_cast<uint32_t>(terms.size())) {
    matrix_.resize(n, Bits(m_ + n));
    coeffs_.resize(m_);
    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j = 0; j < m_; ++j) {
            if (terms[j].first & (1u << i))
                matrix_[i].set(j);
        }
        matrix_[i].set(m_ + i);
    }
    for (uint32_t j = 0; j < m_; ++j)
        coeffs_[j] = terms[j].second;
}

void ParityMatrix::apply_cnot(uint32_t control, uint32_t target) {
    Bits parity_mask(matrix_[0].size());
    for (uint32_t j = 0; j < m_; ++j)
        parity_mask.set(j);
    Bits identity_mask(matrix_[0].size());
    for (uint32_t j = m_; j < m_ + n_; ++j)
        identity_mask.set(j);

    matrix_[control] ^= (matrix_[target] & parity_mask);
    matrix_[target] ^= (matrix_[control] & identity_mask);
}

bool ParityMatrix::is_column_one_hot(uint32_t col) const {
    uint32_t count = 0;
    for (uint32_t i = 0; i < n_; ++i) {
        if (matrix_[i].get(col)) count++;
        if (count > 1) return false;
    }
    return count == 1;
}

uint32_t ParityMatrix::get_one_hot_row(uint32_t col) const {
    for (uint32_t i = 0; i < n_; ++i)
        if (matrix_[i].get(col)) return i;
    return n_;
}

int ParityMatrix::remove_column(uint32_t col) {
    uint32_t col_count = matrix_[0].size();
    for (uint32_t i = 0; i < n_; ++i) {
        Bits new_row(col_count - 1);
        for (uint32_t j = 0; j < col; ++j)
            if (matrix_[i].get(j)) new_row.set(j);
        for (uint32_t j = col + 1; j < col_count; ++j)
            if (matrix_[i].get(j)) new_row.set(j - 1);
        matrix_[i] = new_row;
    }
    int coeff = coeffs_[col];
    coeffs_.erase(coeffs_.begin() + col);
    m_--;
    return coeff;
}

bool ParityMatrix::is_goal_state() const {
    if (m_ != 0) return false;
    for (uint32_t i = 0; i < n_; ++i)
        for (uint32_t j = 0; j < n_; ++j)
            if (matrix_[i].get(j) != (i == j)) return false;
    return true;
}

uint32_t ParityMatrix::count_ones() const {
    uint32_t count = 0;
    for (uint32_t i = 0; i < n_; ++i)
        count += matrix_[i].count();
    return count;
}

uint32_t ParityMatrix::count_ones_in_first_m() const {
    uint32_t count = 0;
    for (uint32_t i = 0; i < n_; ++i)
        for (uint32_t j = 0; j < m_; ++j)
            if (matrix_[i].get(j)) count++;
    return count;
}

uint64_t ParityMatrix::hash() const {
    if (matrix_.empty()) return 0;
    std::hash<Bits> bits_hash;
    uint64_t h = 0;
    for (uint32_t i = 0; i < n_ && i < matrix_.size(); ++i) {
        size_t row_hash = bits_hash(matrix_[i]);
        h ^= static_cast<uint64_t>(row_hash) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
}

ParityMatrix ParityMatrix::copy() const {
    ParityMatrix result(n_, m_);
    result.matrix_ = matrix_;
    result.coeffs_ = coeffs_;
    return result;
}

bool ParityMatrix::operator==(ParityMatrix const& other) const {
    if (n_ != other.n_ || m_ != other.m_) return false;
    if (matrix_ != other.matrix_) return false;
    return coeffs_ == other.coeffs_;
}

bool ParityMatrix::operator<(ParityMatrix const& other) const {
    if (n_ != other.n_) return n_ < other.n_;
    if (m_ != other.m_) return m_ < other.m_;
    for (uint32_t i = 0; i < n_ && i < matrix_.size() && i < other.matrix_.size(); ++i) {
        if (matrix_[i].raw() != other.matrix_[i].raw())
            return matrix_[i].raw() < other.matrix_[i].raw();
    }
    return coeffs_ < other.coeffs_;
}

ParityMatrix ParityMatrix::canonicalize(std::vector<uint32_t>& row_perm,
                                         std::vector<uint32_t>& col_perm) const {
    row_perm.resize(n_);
    for (uint32_t i = 0; i < n_; ++i) row_perm[i] = i;

    std::vector<uint32_t> kept_cols;
    uint32_t col_count = matrix_[0].size();
    for (uint32_t j = 0; j < m_; ++j) {
        bool has_ones = false;
        for (uint32_t i = 0; i < n_; ++i)
            if (matrix_[i].get(j)) { has_ones = true; break; }
        if (has_ones) kept_cols.push_back(j);
    }
    for (uint32_t j = m_; j < col_count; ++j)
        kept_cols.push_back(j);

    uint32_t new_m = static_cast<uint32_t>(kept_cols.size()) - (col_count - m_);
    std::vector<std::pair<std::vector<bool>, uint32_t>> col_sigs;
    for (uint32_t j = 0; j < new_m; ++j) {
        std::vector<bool> col_sig(n_);
        for (uint32_t i = 0; i < n_; ++i)
            col_sig[i] = matrix_[i].get(kept_cols[j]);
        col_sigs.push_back({col_sig, j});
    }
    std::sort(col_sigs.begin(), col_sigs.end());

    ParityMatrix result(n_, new_m);
    result.coeffs_.resize(new_m, 0);
    result.matrix_.resize(n_, Bits(static_cast<uint32_t>(kept_cols.size())));

    col_perm.resize(kept_cols.size());
    for (uint32_t j = 0; j < new_m; ++j) {
        col_perm[j] = kept_cols[col_sigs[j].second];
        for (uint32_t i = 0; i < n_; ++i)
            if (matrix_[i].get(col_perm[j])) result.matrix_[i].set(j);
        if (col_perm[j] < m_)
            result.coeffs_[j] = coeffs_[col_perm[j]];
    }
    for (uint32_t j = new_m; j < kept_cols.size(); ++j) {
        col_perm[j] = kept_cols[j];
        for (uint32_t i = 0; i < n_; ++i)
            if (matrix_[i].get(col_perm[j])) result.matrix_[i].set(j);
    }

    return result;
}

}  // namespace exact_t
