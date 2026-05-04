#include "solver/affine.hpp"
#include "solver/affine-state.hpp"

#include <algorithm>
#include <chrono>
#include <queue>
#include <unordered_set>

namespace exact_t {
namespace affine {

namespace {

bool should_terminate(uint32_t iteration,
                      std::chrono::high_resolution_clock::time_point start,
                      double timeout_s) {
    if (iteration > 100000) return true;
    auto now = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(now - start).count() > timeout_s;
}

uint32_t find_state_index(std::vector<SearchState> const& all_states,
                           SearchState const& target) {
    for (size_t i = 0; i < all_states.size(); ++i)
        if (all_states[i].hash() == target.hash() &&
            all_states[i].g_cost == target.g_cost)
            return static_cast<uint32_t>(i);
    return static_cast<uint32_t>(all_states.size()) - 1;
}

std::vector<uint32_t> compute_active_rows(SearchState const& current, uint32_t n) {
    std::vector<uint32_t> active_rows;
    for (uint32_t i = 0; i < n; ++i) {
        bool active = false;
        for (uint32_t j = 0; j < current.phase_parities.cols() && !active; ++j)
            if (current.phase_parities.get(i, j) &&
                current.phase_parities.count_ones_in_col(j) > 1)
                active = true;
        for (uint32_t j = 0; j < n && !active; ++j)
            if (current.output_parities.get(i, j) != (i == j))
                active = true;
        if (active) active_rows.push_back(i);
    }
    if (active_rows.empty())
        for (uint32_t i = 0; i < n; ++i) active_rows.push_back(i);
    return active_rows;
}

void prune_search_queue(
    std::priority_queue<SearchState, std::vector<SearchState>, SearchStateCompare>& open_set,
    uint32_t max_queue_size) {
    std::vector<SearchState> top;
    while (!open_set.empty() && top.size() < max_queue_size) {
        top.push_back(open_set.top());
        open_set.pop();
    }
    std::priority_queue<SearchState, std::vector<SearchState>, SearchStateCompare> nq;
    for (auto& s : top) nq.push(s);
    open_set = std::move(nq);
}

AffineResult handle_empty_phase(BitMatrix const& target, uint32_t n) {
    AffineResult result;
    BitMatrix current = identity_matrix(n);
    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j = 0; j < n; ++j) {
            if (current.get(i, j) != target.get(i, j)) {
                for (uint32_t k = 0; k < n; ++k) {
                    if (k != i && target.get(k, j)) {
                        result.cnots.push_back({k, i});
                        apply_cnot_to_matrix(current, k, i);
                        break;
                    }
                }
            }
        }
    }
    result.cnot_count = static_cast<uint32_t>(result.cnots.size());
    result.success = true;
    return result;
}

AffineResult handle_initial_goal(Z8Phase const& phase, uint32_t n) {
    AffineResult result;
    result.success = true;
    for (auto const& [parity, coeff] : phase.terms) {
        if (coeff != 0) {
            for (uint32_t q = 0; q < n; ++q) {
                if (parity == (1u << q)) {
                    result.phases.push_back({0, {q, coeff}});
                    break;
                }
            }
        }
    }
    return result;
}

std::vector<std::pair<uint32_t, uint32_t>> reconstruct_cnot_path(
    std::vector<SearchState> const& all_states, SearchState const& goal) {
    uint32_t idx = find_state_index(all_states, goal);
    std::vector<std::pair<uint32_t, uint32_t>> path;
    while (idx < all_states.size() && all_states[idx].cnot.first != UINT32_MAX) {
        path.push_back(all_states[idx].cnot);
        idx = all_states[idx].parent_idx;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

void place_phases_along_path(AffineResult& result, Z8Phase const& phase,
                              uint32_t n,
                              std::vector<std::pair<uint32_t, uint32_t>> const& cnot_path) {
    BitMatrix parity_state = identity_matrix(n);
    std::vector<uint32_t> parities;
    std::vector<int> coeffs;
    for (auto const& [p, c] : phase.terms) {
        if (c != 0) { parities.push_back(p); coeffs.push_back(c); }
    }

    size_t cnot_idx = 0;
    for (size_t term_idx = 0; term_idx < parities.size(); ++term_idx) {
        uint32_t target_parity = parities[term_idx];
        int coeff = coeffs[term_idx];

        while (cnot_idx <= cnot_path.size()) {
            for (uint32_t q = 0; q < n; ++q) {
                uint32_t cp = 0;
                for (uint32_t j = 0; j < n; ++j)
                    if (parity_state.get(q, j)) cp |= (1u << j);
                if (cp == target_parity) {
                    result.phases.push_back({cnot_idx, {q, coeff}});
                    goto next_term;
                }
            }
            if (cnot_idx < cnot_path.size()) {
                auto [ctrl, tgt] = cnot_path[cnot_idx];
                apply_cnot_to_matrix(parity_state, ctrl, tgt);
            }
            ++cnot_idx;
        }
    next_term:;
    }
}

AffineResult astar_search(Z8Phase const& phase, BitMatrix const& target,
                              uint32_t n, AffineParams const& params,
                              double timeout_s) {
    AffineResult result;

    if (phase.terms.empty()) return handle_empty_phase(target, n);

    auto start = std::chrono::high_resolution_clock::now();
    SearchState initial(phase, target, n);

    if (initial.is_goal()) return handle_initial_goal(phase, n);

    std::priority_queue<SearchState, std::vector<SearchState>, SearchStateCompare> open_set;
    std::unordered_set<uint64_t> visited;
    std::vector<SearchState> all_states;

    all_states.push_back(initial);
    open_set.push(initial);
    visited.insert(initial.hash());
    result.states_explored = 1;

    uint32_t iteration = 0;
    while (!open_set.empty()) {
        ++iteration;
        if (should_terminate(iteration, start, timeout_s)) {
            result.states_explored = static_cast<uint32_t>(all_states.size());
            return result;
        }

        SearchState current = open_set.top();
        open_set.pop();
        result.states_explored = static_cast<uint32_t>(all_states.size());

        if (current.is_goal()) {
            auto cnot_path = reconstruct_cnot_path(all_states, current);
            result.cnots = cnot_path;
            result.cnot_count = static_cast<uint32_t>(cnot_path.size());
            result.success = true;
            place_phases_along_path(result, phase, n, cnot_path);
            return result;
        }

        uint32_t parent_idx = find_state_index(all_states, current);
        auto active_rows = compute_active_rows(current, n);

        for (uint32_t ctrl : active_rows) {
            for (uint32_t tgt : active_rows) {
                if (ctrl == tgt) continue;
                SearchState next = current.apply_cnot(ctrl, tgt);
                uint64_t nh = next.hash();
                if (visited.count(nh) == 0) {
                    next.parent_idx = parent_idx;
                    all_states.push_back(next);
                    open_set.push(next);
                    visited.insert(nh);
                    if (visited.size() > params.max_queue_size * 2)
                        prune_search_queue(open_set, params.max_queue_size);
                }
            }
        }
    }
    return result;
}

}  // namespace

AffineResult astar_synthesize(Z8Phase const& phase, BitMatrix const& target,
                                  uint32_t n, AffineParams const& params) {
    return astar_search(phase, target, n, params, params.timeout_s);
}

AffineResult gray_fallback(Z8Phase const& phase, BitMatrix const& target,
                               uint32_t n) {
    AffineResult result;

    std::vector<std::pair<uint32_t, int>> sorted_terms;
    for (auto const& [p, c] : phase.terms)
        if (c != 0) sorted_terms.push_back({p, c});

    std::sort(sorted_terms.begin(), sorted_terms.end(),
              [](auto const& a, auto const& b) {
                  return __builtin_popcount(a.first) < __builtin_popcount(b.first);
              });

    BitMatrix parity_state = identity_matrix(n);

    for (auto const& [target_parity, coeff] : sorted_terms) {
        int pivot = -1;
        for (uint32_t row = 0; row < n; ++row) {
            uint32_t rp = 0;
            for (uint32_t j = 0; j < n; ++j)
                if (parity_state.get(row, j)) rp |= (1u << j);
            if (rp == target_parity) { pivot = static_cast<int>(row); break; }
        }

        if (pivot == -1) {
            for (uint32_t row = 0; row < n; ++row) {
                if (target_parity & (1u << row)) { pivot = static_cast<int>(row); break; }
            }
            if (pivot == -1) continue;

            for (uint32_t row = 0; row < n; ++row) {
                if (row == static_cast<uint32_t>(pivot)) continue;
                uint32_t rp = 0, pp = 0;
                for (uint32_t j = 0; j < n; ++j) {
                    if (parity_state.get(row, j)) rp |= (1u << j);
                    if (parity_state.get(pivot, j)) pp |= (1u << j);
                }
                if ((target_parity & rp) && !(pp & rp)) {
                    result.cnots.push_back({row, static_cast<uint32_t>(pivot)});
                    apply_cnot_to_matrix(parity_state, row, static_cast<uint32_t>(pivot));
                }
            }
        }

        result.phases.push_back({result.cnots.size(), {static_cast<uint32_t>(pivot), coeff}});
    }

    BitMatrix inv_target = invert_matrix(target);
    if (inv_target.rows() == 0) {
        result.cnot_count = static_cast<uint32_t>(result.cnots.size());
        result.success = true;
        return result;
    }

    BitMatrix needed(n, n);
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t j = 0; j < n; ++j) {
            bool val = false;
            for (uint32_t k = 0; k < n; ++k)
                if (inv_target.get(i, k) && parity_state.get(k, j)) val = !val;
            if (val) needed.set(i, j);
        }

    for (uint32_t col = 0; col < n; ++col) {
        uint32_t pivot = n;
        for (uint32_t row = col; row < n; ++row)
            if (needed.get(row, col)) { pivot = row; break; }
        if (pivot == n) continue;

        if (pivot != col) {
            result.cnots.push_back({pivot, col});
            result.cnots.push_back({col, pivot});
            result.cnots.push_back({pivot, col});
            needed.swap_rows(pivot, col);
        }
        for (uint32_t row = 0; row < n; ++row) {
            if (row != col && needed.get(row, col)) {
                result.cnots.push_back({col, row});
                needed.xor_rows(row, col);
            }
        }
    }

    result.cnot_count = static_cast<uint32_t>(result.cnots.size());
    result.success = true;
    return result;
}

Circuit result_to_circuit(AffineResult const& result,
                           Z8Phase const& /*phase*/, uint32_t n) {
    Circuit circuit;
    circuit.set_num_qubits(n);

    auto sorted_phases = result.phases;
    std::sort(sorted_phases.begin(), sorted_phases.end(),
              [](auto const& a, auto const& b) { return a.first < b.first; });

    size_t phase_idx = 0;
    for (size_t k = 0; k <= result.cnots.size(); ++k) {
        while (phase_idx < sorted_phases.size() && sorted_phases[phase_idx].first == k) {
            uint32_t qubit = sorted_phases[phase_idx].second.first;
            int coeff = sorted_phases[phase_idx].second.second;
            circuit.append_phase_gates(qubit, coeff);
            ++phase_idx;
        }
        if (k < result.cnots.size()) {
            auto [ctrl, tgt] = result.cnots[k];
            circuit.cx(ctrl, tgt);
        }
    }

    while (phase_idx < sorted_phases.size()) {
        circuit.append_phase_gates(sorted_phases[phase_idx].second.first,
                                    sorted_phases[phase_idx].second.second);
        ++phase_idx;
    }

    return circuit;
}

}  // namespace affine
}  // namespace exact_t
