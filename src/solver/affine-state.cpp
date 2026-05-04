#include "solver/affine-state.hpp"
#include "core/linear.hpp"

namespace exact_t {
namespace affine {

SearchState::SearchState(Z8Phase const& phase, BitMatrix const& target, uint32_t n)
    : g_cost(0), h_cost(0), parent_idx(UINT32_MAX), cnot({UINT32_MAX, UINT32_MAX}) {
    output_parities = target;

    std::vector<uint32_t> parities;
    for (auto const& [parity, coeff] : phase.terms)
        if (coeff != 0) parities.push_back(parity);

    uint32_t m = static_cast<uint32_t>(parities.size());
    phase_parities = BitMatrix(n, m);

    for (uint32_t j = 0; j < m; ++j) {
        uint32_t parity = parities[j];
        for (uint32_t i = 0; i < n; ++i)
            if (parity & (1u << i)) phase_parities.set(i, j);
    }

    update_heuristic();
}

bool SearchState::is_goal() const {
    for (uint32_t j = 0; j < phase_parities.cols(); ++j)
        if (phase_parities.count_ones_in_col(j) > 1) return false;

    uint32_t n = output_parities.rows();
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t j = 0; j < n; ++j)
            if (output_parities.get(i, j) != (i == j)) return false;

    return true;
}

SearchState SearchState::apply_cnot(uint32_t ctrl, uint32_t tgt) const {
    SearchState next = *this;
    for (uint32_t j = 0; j < next.phase_parities.cols(); ++j)
        if (next.phase_parities.get(ctrl, j))
            next.phase_parities.flip(tgt, j);

    apply_cnot_to_matrix(next.output_parities, ctrl, tgt);

    next.g_cost = g_cost + 1;
    next.cnot = {ctrl, tgt};
    next.update_heuristic();
    return next;
}

uint64_t SearchState::hash() const {
    uint64_t h = 0x9e3779b97f4a7c15ULL;
    for (uint32_t i = 0; i < phase_parities.rows(); ++i)
        for (uint32_t j = 0; j < phase_parities.cols(); ++j)
            if (phase_parities.get(i, j))
                h ^= (static_cast<uint64_t>(i) * 31 + j) + 0x9e3779b9 + (h << 6) + (h >> 2);

    for (uint32_t i = 0; i < output_parities.rows(); ++i)
        for (uint32_t j = 0; j < output_parities.cols(); ++j)
            if (output_parities.get(i, j))
                h ^= (static_cast<uint64_t>(i + 100) * 37 + j) +
                     0x517cc1b727220a95ULL + (h << 6) + (h >> 2);

    return h;
}

uint32_t SearchState::compute_h1() const {
    uint32_t total = 0;
    for (uint32_t j = 0; j < phase_parities.cols(); ++j) {
        uint32_t w = phase_parities.count_ones_in_col(j);
        if (w > 1) total += w - 1;
    }
    return total;
}

uint32_t SearchState::compute_h2() const {
    uint32_t n = output_parities.rows();
    uint32_t cnot_count = 0;
    BitMatrix temp = output_parities;

    for (uint32_t col = 0; col < n; ++col) {
        if (!temp.get(col, col)) {
            for (uint32_t row = col + 1; row < n; ++row) {
                if (temp.get(row, col)) {
                    temp.xor_rows(col, row);
                    ++cnot_count;
                    break;
                }
            }
        }
        if (!temp.get(col, col)) continue;
        for (uint32_t row = 0; row < n; ++row) {
            if (row != col && temp.get(row, col)) {
                temp.xor_rows(row, col);
                ++cnot_count;
            }
        }
    }
    return cnot_count;
}

uint32_t SearchState::compute_heuristic() const {
    return compute_h1() + compute_h2();
}

void SearchState::update_heuristic() {
    h_cost = compute_heuristic();
}

bool SearchState::term_matches_qubit(uint32_t term_idx, uint32_t qubit) const {
    if (term_idx >= phase_parities.cols() || qubit >= phase_parities.rows())
        return false;
    if (!phase_parities.get(qubit, term_idx)) return false;
    return phase_parities.count_ones_in_col(term_idx) == 1;
}

std::vector<uint32_t> SearchState::get_active_columns() const {
    std::vector<uint32_t> active;
    for (uint32_t j = 0; j < phase_parities.cols(); ++j)
        if (phase_parities.count_ones_in_col(j) > 1)
            active.push_back(j);
    return active;
}

}  // namespace affine
}  // namespace exact_t
