/*! \file linear.hpp
 *  \brief Binary matrix operations over GF(2).
 *
 *  Header-only: BitMatrix, Gaussian elimination, kernel, solve.
 */

#pragma once

#include "core/bits.hpp"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

namespace exact_t {

/*! \brief Binary matrix with row operations over GF(2). */
class BitMatrix {
  public:
    BitMatrix() = default;

    /*! \brief Construct zero matrix of given dimensions. */
    BitMatrix(uint32_t rows, uint32_t cols) : rows_(rows), cols_(cols) {
        data_.reserve(rows);
        for (uint32_t i = 0; i < rows; ++i)
            data_.emplace_back(cols);
    }

    uint32_t rows() const { return rows_; }
    uint32_t cols() const { return cols_; }
    std::vector<Bits> const& data() const { return data_; }

    Bits& operator[](uint32_t row) { return data_[row]; }
    Bits const& operator[](uint32_t row) const { return data_[row]; }

    bool get(uint32_t row, uint32_t col) const { return data_[row].get(col); }
    void set(uint32_t row, uint32_t col) { data_[row].set(col); }
    void flip(uint32_t row, uint32_t col) { data_[row].flip(col); }

    /*! \brief Swap two rows. */
    void swap_rows(uint32_t r1, uint32_t r2) {
        if (r1 != r2) std::swap(data_[r1], data_[r2]);
    }

    /*! \brief XOR src row into dst row. */
    void xor_rows(uint32_t dst, uint32_t src) { data_[dst] ^= data_[src]; }

    /*! \brief Return transposed matrix. */
    BitMatrix transpose() const {
        BitMatrix result(cols_, rows_);
        for (uint32_t i = 0; i < rows_; ++i)
            for (uint32_t j = 0; j < cols_; ++j)
                if (get(i, j)) result.set(j, i);
        return result;
    }

    /*! \brief Count total ones in matrix. */
    uint32_t count_ones() const {
        uint32_t count = 0;
        for (auto const& row : data_) count += row.count();
        return count;
    }

    /*! \brief Count ones in column. */
    uint32_t count_ones_in_col(uint32_t col) const {
        uint32_t count = 0;
        for (uint32_t i = 0; i < rows_; ++i)
            if (get(i, col)) count++;
        return count;
    }

    /*! \brief Check if column has exactly one 1. */
    bool is_column_one_hot(uint32_t col) const {
        return count_ones_in_col(col) == 1;
    }

    /*! \brief Find row index with 1 in one-hot column. */
    uint32_t get_one_hot_row(uint32_t col) const {
        for (uint32_t i = 0; i < rows_; ++i)
            if (get(i, col)) return i;
        return rows_;
    }

  private:
    uint32_t rows_{0};
    uint32_t cols_{0};
    std::vector<Bits> data_;
};

/*! \brief Result of Gaussian elimination. */
struct GaussianResult {
    BitMatrix reduced;
    std::vector<int> pivot_cols;
    uint32_t rank{0};
    BitMatrix transform;
};

/*! \brief Perform Gaussian elimination to row echelon form. */
inline GaussianResult gaussian_eliminate(BitMatrix const& mat,
                                          bool compute_transform = false) {
    GaussianResult result;
    result.reduced = mat;
    result.rank = 0;
    result.pivot_cols.resize(mat.rows(), -1);

    if (compute_transform) {
        result.transform = BitMatrix(mat.rows(), mat.rows());
        for (uint32_t i = 0; i < mat.rows(); ++i)
            result.transform.set(i, i);
    }

    uint32_t pivot_row = 0;
    for (uint32_t col = 0; col < mat.cols() && pivot_row < mat.rows(); ++col) {
        uint32_t found_row = pivot_row;
        while (found_row < mat.rows() && !result.reduced.get(found_row, col))
            found_row++;
        if (found_row == mat.rows()) continue;

        if (found_row != pivot_row) {
            result.reduced.swap_rows(pivot_row, found_row);
            if (compute_transform)
                result.transform.swap_rows(pivot_row, found_row);
        }

        result.pivot_cols[pivot_row] = static_cast<int>(col);

        for (uint32_t r = 0; r < mat.rows(); ++r) {
            if (r != pivot_row && result.reduced.get(r, col)) {
                result.reduced.xor_rows(r, pivot_row);
                if (compute_transform)
                    result.transform.xor_rows(r, pivot_row);
            }
        }
        pivot_row++;
    }

    result.rank = pivot_row;
    return result;
}

/*! \brief Find basis vectors for matrix kernel (null space). */
inline std::vector<Bits> find_kernel_basis(BitMatrix const& mat) {
    GaussianResult ge = gaussian_eliminate(mat);
    std::vector<Bits> basis;
    Bits is_pivot(mat.cols());
    for (uint32_t r = 0; r < ge.rank; ++r)
        if (ge.pivot_cols[r] >= 0) is_pivot.set(ge.pivot_cols[r]);

    for (uint32_t c = 0; c < mat.cols(); ++c) {
        if (!is_pivot.get(c)) {
            Bits kernel_vec(mat.cols());
            kernel_vec.set(c);
            for (uint32_t r = 0; r < ge.rank; ++r)
                if (ge.pivot_cols[r] >= 0 && ge.reduced.get(r, c))
                    kernel_vec.set(ge.pivot_cols[r]);
            basis.push_back(kernel_vec);
        }
    }
    return basis;
}

/*! \brief Compute matrix rank. */
inline uint32_t compute_rank(BitMatrix const& mat) {
    return gaussian_eliminate(mat).rank;
}

/*! \brief Create n x n identity matrix. */
inline BitMatrix identity_matrix(uint32_t n) {
    BitMatrix mat(n, n);
    for (uint32_t i = 0; i < n; ++i) mat.set(i, i);
    return mat;
}

/*! \brief Apply CNOT transformation to matrix rows. */
inline void apply_cnot_to_matrix(BitMatrix& mat, uint32_t control,
                                  uint32_t target) {
    mat.xor_rows(target, control);
}

/*! \brief Multiply two binary matrices over GF(2). */
inline BitMatrix multiply_matrices(BitMatrix const& A, BitMatrix const& B) {
    assert(A.cols() == B.rows());
    BitMatrix C(A.rows(), B.cols());
    for (uint32_t i = 0; i < A.rows(); ++i)
        for (uint32_t j = 0; j < B.cols(); ++j) {
            bool sum = false;
            for (uint32_t k = 0; k < A.cols(); ++k)
                sum ^= (A.get(i, k) && B.get(k, j));
            if (sum) C.set(i, j);
        }
    return C;
}

/*! \brief Compute matrix inverse (empty if singular). */
inline BitMatrix invert_matrix(BitMatrix const& mat) {
    assert(mat.rows() == mat.cols());
    uint32_t n = mat.rows();
    BitMatrix augmented(n, 2 * n);
    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j = 0; j < n; ++j)
            if (mat.get(i, j)) augmented.set(i, j);
        augmented.set(i, n + i);
    }

    for (uint32_t col = 0; col < n; ++col) {
        uint32_t pivot = col;
        while (pivot < n && !augmented.get(pivot, col)) pivot++;
        if (pivot == n) return BitMatrix();
        augmented.swap_rows(col, pivot);
        for (uint32_t r = 0; r < n; ++r)
            if (r != col && augmented.get(r, col))
                augmented.xor_rows(r, col);
    }

    BitMatrix inverse(n, n);
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t j = 0; j < n; ++j)
            if (augmented.get(i, n + j)) inverse.set(i, j);
    return inverse;
}

}  // namespace exact_t
