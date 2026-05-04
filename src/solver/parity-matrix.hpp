/*! \file parity-matrix.hpp
 *  \brief Matrix representation for parity network synthesis.
 */

#pragma once

#include "core/bits.hpp"
#include <cstdint>
#include <vector>

namespace exact_t {

/*! \brief Matrix for A* parity network synthesis.
 *
 *  Represents n qubits x (m parities + n identity) columns.
 *  Used in A* search for optimal CNOT synthesis.
 */
class ParityMatrix {
  public:
    ParityMatrix() = default;
    ParityMatrix(uint32_t n, uint32_t m);
    ParityMatrix(std::vector<std::pair<uint32_t, int>> const& terms, uint32_t n);

    uint32_t n() const { return n_; }
    uint32_t m() const { return m_; }
    std::vector<Bits> const& matrix() const { return matrix_; }
    std::vector<int> const& coeffs() const { return coeffs_; }
    Bits& row(uint32_t i) { return matrix_[i]; }
    Bits const& row(uint32_t i) const { return matrix_[i]; }
    int& coeff(uint32_t col) { return coeffs_[col]; }

    /*! \brief Apply CNOT to matrix rows. */
    void apply_cnot(uint32_t control, uint32_t target);

    /*! \brief Check if column has exactly one 1. */
    bool is_column_one_hot(uint32_t col) const;

    /*! \brief Get row index with 1 in one-hot column. */
    uint32_t get_one_hot_row(uint32_t col) const;

    /*! \brief Remove column (parity synthesized). Returns coefficient. */
    int remove_column(uint32_t col);

    /*! \brief Check if all parities synthesized and identity restored. */
    bool is_goal_state() const;

    /*! \brief Count ones in entire matrix. */
    uint32_t count_ones() const;

    /*! \brief Count ones in parity columns only. */
    uint32_t count_ones_in_first_m() const;

    /*! \brief Compute hash for deduplication. */
    uint64_t hash() const;

    /*! \brief Create deep copy. */
    ParityMatrix copy() const;
    bool operator==(ParityMatrix const& other) const;
    bool operator<(ParityMatrix const& other) const;

    /*! \brief Compute canonical form. */
    ParityMatrix canonicalize(std::vector<uint32_t>& row_perm,
                              std::vector<uint32_t>& col_perm) const;

  private:
    uint32_t n_{0};
    uint32_t m_{0};
    std::vector<Bits> matrix_;
    std::vector<int> coeffs_;
};

/*! \brief Hash functor for ParityMatrix. */
struct ParityMatrixHash {
    std::size_t operator()(ParityMatrix const& pmat) const {
        return pmat.hash();
    }
};

}  // namespace exact_t
