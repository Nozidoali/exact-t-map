#include "solver/parity-solver.hpp"

#include <algorithm>
#include <chrono>

namespace exact_t {

namespace {

double elapsed_ms(std::chrono::high_resolution_clock::time_point start) {
    auto now = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(now - start).count();
}

}  // namespace

ParitySolver::ParitySolver(uint32_t n, ParitySolverParam const& params)
    : n_(n), params_(params) {
    if (!params_.cache_path.empty())
        graph_.load_file(params_.cache_path);
}

Circuit ParitySolver::gray_synth(ParityMatrix const& pmat) const {
    auto start = std::chrono::high_resolution_clock::now();

    Circuit c;
    c.set_num_qubits(pmat.n());
    ParityMatrix pm = pmat.copy();

    for (uint32_t col = 0; col < pm.m(); ++col) {
        int pivot = -1;
        for (uint32_t row = 0; row < pm.n(); ++row)
            if (pm.row(row).get(col)) { pivot = static_cast<int>(row); break; }
        if (pivot == -1) continue;

        for (uint32_t row = 0; row < pm.n(); ++row) {
            if (row != static_cast<uint32_t>(pivot) && pm.row(row).get(col)) {
                c.cx(row, pivot);
                pm.apply_cnot(row, pivot);
            }
        }
        if (pm.is_column_one_hot(col)) {
            uint32_t r = pm.get_one_hot_row(col);
            c.append_phase_gates(r, pm.coeffs()[col]);
        }
    }

    for (uint32_t col = 0; col < pm.n(); ++col) {
        int pivot = -1;
        for (uint32_t row = col; row < pm.n(); ++row)
            if (pm.row(row).get(pm.m() + col)) { pivot = static_cast<int>(row); break; }
        if (pivot == -1) continue;
        if (static_cast<uint32_t>(pivot) != col) {
            c.cx(pivot, col); pm.apply_cnot(pivot, col);
            c.cx(col, pivot); pm.apply_cnot(col, pivot);
            c.cx(pivot, col); pm.apply_cnot(pivot, col);
        }
        for (uint32_t row = 0; row < pm.n(); ++row) {
            if (row != col && pm.row(row).get(pm.m() + col)) {
                c.cx(col, row);
                pm.apply_cnot(col, row);
            }
        }
    }

    stats_.gray_synth_time_ms += elapsed_ms(start);
    return c;
}

Circuit ParitySolver::synthesize(Z8Phase const& z8_phase) const {
    stats_.total_calls++;
    auto start = std::chrono::high_resolution_clock::now();

    if (z8_phase.empty()) {
        Circuit c;
        c.set_num_qubits(z8_phase.n);
        return c;
    }

    auto sorted_terms = z8_phase.to_vector();
    ParityMatrix pmat(sorted_terms, z8_phase.n);
    Circuit result = astar_synth(pmat);

    stats_.total_time_ms += elapsed_ms(start);
    return result;
}

Circuit ParitySolver::astar_synth(ParityMatrix const& initial_pmat) const {
    auto start = std::chrono::high_resolution_clock::now();

    AStarContext ctx;

    PSearchState initial;
    initial.pmat = initial_pmat;
    initial.cost = 0;
    initial.heuristic = initial.pmat.count_ones();
    ctx.pq.push(std::move(initial));

    uint32_t iterations = 0;
    Circuit found_result;
    bool found = false;

    while (!ctx.pq.empty() && iterations < params_.max_iterations) {
        iterations++;
        PSearchState current = ctx.pq.top();
        ctx.pq.pop();

        std::vector<uint32_t> rp, cp;
        ParityMatrix current_canonical = current.pmat.canonicalize(rp, cp);
        if (ctx.visited.count(current_canonical)) continue;
        ctx.visited.insert(current_canonical);

        ParityMatrix canon_copy = current_canonical;
        if (process_phases(current.pmat, canon_copy, ctx)) {
            PSearchState next;
            next.pmat = current.pmat;
            next.cost = current.cost;
            next.heuristic = current.pmat.count_ones();
            ctx.pq.push(std::move(next));
            continue;
        }

        if (current.pmat.is_goal_state()) {
            found_result = build_path(current_canonical, initial_pmat, ctx);
            found = true;
            break;
        }

        Graph& graph = graph_.get_graph(initial_pmat.n());
        graph_.record_lookup();
        stats_.cache_lookups++;

        uint32_t graph_node_id = graph.find_node(current_canonical);
        if (graph_node_id != UINT32_MAX) {
            GraphNode const* gn = graph.get_node(graph_node_id);
            if (gn && gn->edge.next_id != UINT32_MAX) {
                graph_.record_hit();
                stats_.cache_hits++;
                found_result = build_path(current_canonical, initial_pmat, ctx);
                Circuit gc = graph.extract_circuit(graph_node_id, current.pmat);
                std::vector<uint32_t> identity(gc.num_qubits());
    for (uint32_t i = 0; i < gc.num_qubits(); ++i) identity[i] = i;
    found_result.append(gc, identity);
                found = true;
                break;
            }
            graph_.record_miss();
            stats_.cache_misses++;
        } else {
            graph_.record_miss();
            stats_.cache_misses++;
        }

        auto candidates = gen_cnot_candidates(current.pmat, initial_pmat.n());
        uint32_t n_explore = std::min(params_.n_branches,
                                       static_cast<uint32_t>(candidates.size()));
        for (uint32_t i = 0; i < n_explore; ++i) {
            auto [ctrl, tgt, score] = candidates[i];
            if (score < 0) continue;
            expand_cnot(current.pmat, current_canonical, current.cost, ctrl, tgt, ctx);
        }
    }

    stats_.astar_time_ms += elapsed_ms(start);

    if (!found)
        return gray_synth(initial_pmat);
    return found_result;
}

bool ParitySolver::process_phases(ParityMatrix& pmat,
                                   ParityMatrix& canon,
                                   AStarContext& ctx) const {
    for (uint32_t col = 0; col < pmat.m(); ++col) {
        if (pmat.is_column_one_hot(col)) {
            uint32_t row = pmat.get_one_hot_row(col);
            uint32_t col_idx = col;
            pmat.remove_column(col);

            std::vector<uint32_t> rp2, cp2;
            ParityMatrix next_canon = pmat.canonicalize(rp2, cp2);

            if (!ctx.parent_map.count(next_canon)) {
                ParentInfo info;
                info.parent_canonical = canon;
                info.is_phase_gate = true;
                info.phase_row = row;
                info.phase_col = col_idx;
                info.row_perm = rp2;
                info.col_perm = cp2;
                ctx.parent_map[next_canon] = info;
            }
            canon = next_canon;
            return true;
        }
    }
    return false;
}

Circuit ParitySolver::build_path(ParityMatrix const& goal_canon,
                                  ParityMatrix const& initial_pmat,
                                  AStarContext const& ctx) const {
    Graph& graph = graph_.get_graph(initial_pmat.n());

    std::vector<ParentInfo const*> path;
    ParityMatrix state = goal_canon;
    uint32_t goal_id = graph.add_node(goal_canon);
    graph.set_node_cost(goal_id, 0);

    int steps = 0;
    while (ctx.parent_map.count(state)) {
        ParentInfo const& info = ctx.parent_map.at(state);
        path.push_back(&info);
        uint32_t curr_id = graph.add_node(state);
        uint32_t parent_id = graph.add_node(info.parent_canonical);
        graph.set_edge(parent_id, curr_id, info.is_phase_gate,
                       info.control, info.target,
                       info.phase_row, info.phase_col,
                       info.row_perm, info.col_perm);
        graph.set_node_cost(curr_id, steps);
        graph.set_node_cost(parent_id, steps + 1);
        steps++;
        state = info.parent_canonical;
    }

    std::reverse(path.begin(), path.end());

    Circuit circuit;
    circuit.set_num_qubits(initial_pmat.n());
    ParityMatrix pmat = initial_pmat;
    for (auto const* info : path) {
        if (info->is_phase_gate) {
            if (info->phase_col < pmat.m() &&
                pmat.is_column_one_hot(info->phase_col) &&
                pmat.get_one_hot_row(info->phase_col) == info->phase_row) {
                int coeff = pmat.remove_column(info->phase_col);
                circuit.append_phase_gates(info->phase_row, coeff);
            }
        } else {
            if (info->control != info->target &&
                info->control < pmat.n() && info->target < pmat.n()) {
                circuit.cx(info->control, info->target);
                pmat.apply_cnot(info->control, info->target);
            }
        }
    }
    return circuit;
}

std::vector<std::tuple<uint32_t, uint32_t, int>>
ParitySolver::gen_cnot_candidates(ParityMatrix const& pmat, uint32_t n) const {
    std::vector<std::tuple<uint32_t, uint32_t, int>> candidates;
    for (uint32_t ctrl = 0; ctrl < n; ++ctrl) {
        for (uint32_t tgt = 0; tgt < n; ++tgt) {
            if (ctrl == tgt) continue;
            ParityMatrix test = pmat.copy();
            uint32_t ones_before = test.count_ones();
            test.apply_cnot(ctrl, tgt);
            int score = static_cast<int>(ones_before) - static_cast<int>(test.count_ones());
            candidates.emplace_back(ctrl, tgt, score);
        }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](auto const& a, auto const& b) {
                  return std::get<2>(a) > std::get<2>(b);
              });
    return candidates;
}

void ParitySolver::expand_cnot(ParityMatrix const& current_pmat,
                                ParityMatrix const& current_canonical,
                                uint32_t cost,
                                uint32_t ctrl, uint32_t tgt,
                                AStarContext& ctx) const {
    PSearchState next;
    next.pmat = current_pmat.copy();
    next.pmat.apply_cnot(ctrl, tgt);
    next.cost = cost + 1;
    next.heuristic = next.pmat.count_ones();

    std::vector<uint32_t> rp2, cp2;
    ParityMatrix next_canon = next.pmat.canonicalize(rp2, cp2);

    if (ctx.visited.count(next_canon)) return;

    if (!ctx.parent_map.count(next_canon)) {
        ParentInfo info;
        info.parent_canonical = current_canonical;
        info.is_phase_gate = false;
        info.control = ctrl;
        info.target = tgt;
        info.dist_to_goal = next.heuristic + next.cost;
        info.row_perm = rp2;
        info.col_perm = cp2;
        ctx.parent_map[next_canon] = info;
    }

    ctx.pq.push(std::move(next));
}

void ParitySolver::load_cache(std::string const& path) {
    if (!path.empty()) graph_.load_file(path);
}

void ParitySolver::save_cache(std::string const& path) const {
    if (!path.empty()) graph_.save_file(path);
}

}  // namespace exact_t
