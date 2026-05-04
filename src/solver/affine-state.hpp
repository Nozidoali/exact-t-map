/*! \file affine-state.hpp
 *  \brief Search state for Affine A* optimization.
 */

#pragma once

#include "core/linear.hpp"
#include "solver/solver.hpp"

#include <cstdint>
#include <utility>
#include <vector>

namespace exact_t {
namespace affine {

/*! \brief Search state for A* phase polynomial optimization.
 *
 *  Represents combined [P | O] matrix where:
 *  - P: n x m phase parity matrix (columns = remaining phase terms)
 *  - O: n x n output linear transformation
 */
struct SearchState {
    BitMatrix phase_parities;
    BitMatrix output_parities;
    uint32_t g_cost{0};
    uint32_t h_cost{0};
    uint32_t parent_idx{UINT32_MAX};
    std::pair<uint32_t, uint32_t> cnot{UINT32_MAX, UINT32_MAX};

    SearchState() = default;

    /*! \brief Construct initial state from phase polynomial and target. */
    SearchState(Z8Phase const& phase, BitMatrix const& target, uint32_t n);

    /*! \brief Total estimated cost f = g + h. */
    uint32_t f_cost() const { return g_cost + h_cost; }

    /*! \brief Check if goal: single-qubit phases and identity output. */
    bool is_goal() const;

    /*! \brief Apply CNOT to create new state. */
    SearchState apply_cnot(uint32_t ctrl, uint32_t tgt) const;

    /*! \brief Compute hash for visited set. */
    uint64_t hash() const;

    /*! \brief Hamming weight heuristic for phase matrix. */
    uint32_t compute_h1() const;

    /*! \brief Gaussian elimination heuristic for output. */
    uint32_t compute_h2() const;

    /*! \brief Combined heuristic. */
    uint32_t compute_heuristic() const;

    /*! \brief Update h_cost. */
    void update_heuristic();

    uint32_t num_qubits() const { return output_parities.rows(); }
    uint32_t num_terms() const { return phase_parities.cols(); }

    /*! \brief Check if phase term matches single qubit. */
    bool term_matches_qubit(uint32_t term_idx, uint32_t qubit) const;

    /*! \brief Get columns with weight > 1. */
    std::vector<uint32_t> get_active_columns() const;
};

/*! \brief Comparison for priority queue (min-heap by f-cost). */
struct SearchStateCompare {
    bool operator()(SearchState const& a, SearchState const& b) const {
        if (a.f_cost() != b.f_cost()) return a.f_cost() > b.f_cost();
        return a.g_cost < b.g_cost;
    }
};

}  // namespace affine
}  // namespace exact_t
