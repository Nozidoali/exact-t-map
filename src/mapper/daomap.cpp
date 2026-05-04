#include "mapper/daomap.hpp"

#include <algorithm>
#include <chrono>
#include <limits>
#include <queue>
#include <set>
#include <sstream>
#include <unordered_set>

namespace exact_t {

DaoMapOptimizer::DaoMapOptimizer(xag_network const& xag, topo_view& topo,
                                   CutManager& cut_manager,
                                   CliffordTSolver& solver,
                                   Params const& params)
    : xag_(xag), topo_(topo), cut_manager_(cut_manager), solver_(solver),
      params_(params),
      cost_evaluator_(solver, params.max_solver_n, params.use_clean_ancilla,
                      params.toffoli_mapping),
      selected_match_(xag_), AF_(xag_), ref_(xag_),
      est_ref_(xag_), occ_(xag_), need_(xag_),
      prev_selected_(xag_),
      sink_(make_trace_sink(params.trace_out_path)) {}

void DaoMapOptimizer::emit_phase_marker(std::string const& phase,
                                         std::string const& event,
                                         double time_ms) {
    std::ostringstream oss;
    oss << "{\"_kind\":\"phase_marker\",\"phase\":\"" << phase
        << "\",\"event\":\"" << event << "\""
        << ",\"time_ms\":" << time_ms << "}";
    sink_->emit(oss.str());
}

void DaoMapOptimizer::emit_round_events(std::string const& phase, uint32_t round) {
    uint32_t nodes_changed = 0;

    xag_.foreach_node([&](auto n) {
        if (xag_.is_constant(n) || xag_.is_pi(n)) return;
        CutMatch const* m = selected(n);
        if (m == nullptr) return;

        bool cut_changed = have_prev_selected_
            ? (prev_selected_[n] != selected_match_[n])
            : false;
        if (cut_changed) ++nodes_changed;

        std::ostringstream oss;
        oss << "{\"_kind\":\"node\",\"phase\":\"" << phase << "\",\"round\":" << round
            << ",\"node\":" << n
            << ",\"selected_cut_idx\":" << selected_match_[n]
            << ",\"children\":[";
        for (size_t i = 0; i < m->children.size(); ++i) {
            if (i) oss << ',';
            oss << m->children[i];
        }
        oss << "],\"outputs\":[";
        for (size_t i = 0; i < m->outputs.size(); ++i) {
            if (i) oss << ',';
            oss << m->outputs[i];
        }
        oss << "],\"cut_size\":" << m->children.size()
            << ",\"num_outputs\":" << m->outputs.size()
            << ",\"degree\":" << m->degree()
            << ",\"is_self_cut\":" << (m->is_self_cut() ? "true" : "false")
            << ",\"cost\":{\"t\":" << m->c.t_count
            << ",\"cnot\":" << m->c.cnot_count
            << ",\"q\":" << m->c.q_count
            << ",\"s\":" << m->c.s_count
            << ",\"z\":" << m->c.z_count << "}"
            << ",\"cut_changed\":" << (cut_changed ? "true" : "false")
            << "}";
        sink_->emit(oss.str());
    });

    Cost total_cost = compute_mapped_total_cost();

    std::ostringstream summary;
    summary << "{\"_kind\":\"round_summary\",\"phase\":\"" << phase
            << "\",\"round\":" << round
            << ",\"total_t\":" << total_cost.t_count
            << ",\"total_cnot\":" << total_cost.cnot_count
            << ",\"nodes_changed_since_prev\":" << nodes_changed << "}";
    sink_->emit(summary.str());

    xag_.foreach_node([&](auto n) { prev_selected_[n] = selected_match_[n]; });
    have_prev_selected_ = true;
}

void DaoMapOptimizer::optimize() {
    reset_state();
    topo_.foreach_node([&](auto n) {
        AF_[n] = xag_.is_constant(n) || xag_.is_pi(n) ? 0.0 : params_.area_flow_init;
    });

    auto t0 = std::chrono::high_resolution_clock::now();
    cut_manager_.cut_enum_bilinear();
    double cut_enum_ms = std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - t0).count();
    {
        std::ostringstream oss;
        oss << "{\"_kind\":\"phase_marker\",\"phase\":\"cut_enum\""
            << ",\"event\":\"done\""
            << ",\"total_cuts\":" << cut_manager_.get_stats().total_cuts
            << ",\"time_ms\":" << cut_enum_ms << "}";
        sink_->emit(oss.str());
    }

    topo_.foreach_node([&](auto n) {
        if (xag_.is_constant(n) || xag_.is_pi(n)) return;
        for (uint32_t ci : cut_manager_.get_candidate_matches()[n])
            cut_manager_.get_cut_matches()[ci].c =
                cost_evaluator_.evaluate(cut_manager_.get_cut_matches()[ci]);
    });

    initialize_po_references();
    reverse_topo_order_ = compute_reverse_topological_order();

    std::unordered_map<xag_node, uint32_t> fanout_counts;
    topo_.foreach_node([&](auto n) {
        if (xag_.is_constant(n) || xag_.is_pi(n)) {
            need_[n] = 0;
        } else {
            xag_.foreach_fanin(n, [&](auto fi) {
                fanout_counts[xag_.get_node(fi)]++;
            });
        }
    });
    topo_.foreach_node([&](auto n) {
        if (!xag_.is_constant(n) && !xag_.is_pi(n))
            need_[n] = fanout_counts.count(n) ? fanout_counts[n] : 0;
    });

    for (uint32_t it = 0; it < params_.area_flow_rounds; ++it) {
        if (it == 0) initialize_area_flow();
        else compute_area_flow();
        evaluate_and_select_cuts(it);
        emit_round_events("area_flow", it);
        if (it < params_.area_flow_rounds - 1) {
            update_references_from_exact_mapping(it);
            compute_area_flow();
        }
    }

    if (params_.exact_rounds > 0) {
        topo_.foreach_node([&](auto n) { ref_[n] = 0; });
        xag_.foreach_po([&](auto po) {
            xag_node po_node = xag_.get_node(po);
            if (selected(po_node))
                cut_update_refs(selected_match_[po_node], +1);
        });

        for (uint32_t r = 0; r < params_.exact_rounds; ++r) {
            evaluate_and_select_cuts_exact();
            emit_round_events("exact", r);
            auto exact_refs = compute_mapped_refs();
            double coef = 1.0 / ((params_.area_flow_rounds + r + 2.0) *
                                  (params_.area_flow_rounds + r + 2.0));
            smooth_est_refs(coef, exact_refs);
            compute_area_flow();
        }
    }
}

CutMatch const* DaoMapOptimizer::selected(xag_node n) const {
    uint32_t ci = selected_match_[n];
    if (ci == 0 || ci >= cut_manager_.get_cut_matches().size()) return nullptr;
    return &cut_manager_.get_cut_matches()[ci];
}

Cost DaoMapOptimizer::compute_mapped_total_cost() const {
    std::queue<xag_node> q;
    std::unordered_set<xag_node> visited;
    Cost total = {0, 0, 0, 0, 0};

    xag_.foreach_po([&](auto po) { q.push(xag_.get_node(po)); });

    while (!q.empty()) {
        xag_node n = q.front(); q.pop();
        if (visited.count(n)) continue;
        visited.insert(n);
        if (xag_.is_constant(n) || xag_.is_pi(n)) continue;

        CutMatch const* match = selected(n);
        if (match) {
            total = total + match->c;
            if (match->outputs.size() > 1)
                for (xag_node o : match->outputs) visited.insert(o);
            for (auto const& child : match->children) q.push(child);
        }
    }
    return total;
}

void DaoMapOptimizer::reset_state() {
    topo_.foreach_node([&](auto n) {
        ref_[n] = 0; est_ref_[n] = 1; occ_[n] = 0; need_[n] = 0; AF_[n] = 0.0;
    });
}

void DaoMapOptimizer::initialize_po_references() {
    xag_.foreach_po([&](auto po) { ref_[xag_.get_node(po)]++; });
}

std::vector<xag_node> DaoMapOptimizer::compute_reverse_topological_order() const {
    std::vector<xag_node> rto;
    topo_.foreach_node([&](auto n) {
        if (!xag_.is_constant(n) && !xag_.is_pi(n)) rto.push_back(n);
    });
    std::reverse(rto.begin(), rto.end());
    return rto;
}

void DaoMapOptimizer::initialize_area_flow() {
    topo_.foreach_node([&](auto n) {
        if (xag_.is_constant(n) || xag_.is_pi(n)) { AF_[n] = 0.0; return; }
        double ac = params_.area_flow_init, ffs = 0.0;
        xag_.foreach_fanin(n, [&](auto fi) { ffs += AF_[xag_.get_node(fi)]; });
        AF_[n] = (ac + ffs) / compute_refs_eff(ref_[n], need_[n], params_.alpha);
    });
}

void DaoMapOptimizer::compute_area_flow() {
    topo_.foreach_node([&](auto n) {
        if (xag_.is_constant(n) || xag_.is_pi(n)) { AF_[n] = 0.0; return; }
        double ac = 0.0, ffs = 0.0;
        CutMatch const* match = selected(n);
        if (match) {
            uint32_t no = match->outputs.empty() ? 1u : static_cast<uint32_t>(match->outputs.size());
            ac = static_cast<double>(match->c.t_count) / static_cast<double>(no);
            for (auto const& child : match->children) ffs += AF_[child];
        } else {
            ac = params_.area_flow_init;
            xag_.foreach_fanin(n, [&](auto fi) { ffs += AF_[xag_.get_node(fi)]; });
        }
        AF_[n] = (ac + ffs) / compute_refs_eff(ref_[n], need_[n], params_.alpha);
    });
}

void DaoMapOptimizer::smooth_est_refs(
    double coef,
    std::unordered_map<xag_node, uint32_t> const& exact_refs) {
    topo_.foreach_node([&](auto n) {
        uint32_t er = exact_refs.count(n) ? exact_refs.at(n) : 0;
        double ne = coef * static_cast<double>(est_ref_[n]) +
                    (1.0 - coef) * static_cast<double>(er);
        est_ref_[n] = std::max(1u, static_cast<uint32_t>(ne));
    });
}

void DaoMapOptimizer::update_references_from_exact_mapping(uint32_t iteration) {
    auto exact_refs = compute_mapped_refs();
    double coef = 1.0 / ((iteration + 2.0) * (iteration + 2.0));
    smooth_est_refs(coef, exact_refs);
    topo_.foreach_node([&](auto n) { ref_[n] = 0; occ_[n] = 0; need_[n] = 0; });
    initialize_po_references();
}

void DaoMapOptimizer::evaluate_and_select_cuts(uint32_t) {
    topo_.foreach_node([&](auto n) { occ_[n] = 0; });
    std::unordered_map<xag_node, uint32_t> current_occ;
    covered_.clear();

    for (auto n : reverse_topo_order_) {
        if (covered_.count(n)) continue;
        if (cut_manager_.get_candidate_matches()[n].empty()) continue;
        double best_t = std::numeric_limits<double>::infinity();
        uint32_t best_cx = UINT32_MAX, best_q = UINT32_MAX, best_idx = 0;
        bool has_valid = false;

        for (uint32_t ci : cut_manager_.get_candidate_matches()[n]) {
            auto const& match = cut_manager_.get_cut_matches()[ci];
            if (match.is_self_cut()) continue;

            double ac = static_cast<double>(match.c.t_count);
            double ffs = 0.0, olp = 0.0;
            for (auto& u : match.children) {
                ffs += AF_[u] / compute_refs_eff(ref_[u], need_[u], params_.alpha);
                uint32_t oc = current_occ.count(u) ? current_occ.at(u) : 0;
                olp += params_.beta * std::max(0.0, static_cast<double>(oc) - 1.0);
            }
            double total = ac + ffs + olp;

            if (total < best_t || (total == best_t && match.c.cnot_count < best_cx) ||
                (total == best_t && match.c.cnot_count == best_cx && match.c.q_count < best_q)) {
                best_t = total; best_cx = match.c.cnot_count;
                best_q = match.c.q_count; best_idx = ci; has_valid = true;
            }
        }

        if (has_valid) {
            selected_match_[n] = best_idx;
            auto const& best_match = cut_manager_.get_cut_matches()[best_idx];
            if (best_match.outputs.size() > 1) {
                for (xag_node o : best_match.outputs) {
                    if (o != n) {
                        selected_match_[o] = best_idx;
                        covered_.insert(o);
                        for (auto& u : best_match.children)
                            ref_[u] += 1;
                    }
                }
            }
            for (auto& u : best_match.children) {
                ref_[u] += 1; current_occ[u]++; occ_[u] = current_occ[u];
            }
        }
    }
}

void DaoMapOptimizer::evaluate_and_select_cuts_exact() {
    covered_.clear();

    for (auto n : reverse_topo_order_) {
        if (covered_.count(n)) continue;
        if (cut_manager_.get_candidate_matches()[n].empty()) continue;

        uint32_t best_ec = UINT32_MAX, best_idx = 0;
        bool has_valid = false;
        uint32_t old_idx = selected_match_[n];

        if (selected(n) && ref_[n] > 0)
            cut_update_refs(old_idx, -1);

        for (uint32_t ci : cut_manager_.get_candidate_matches()[n]) {
            auto const& match = cut_manager_.get_cut_matches()[ci];
            if (match.is_self_cut()) continue;
            uint32_t ec = cut_update_refs(ci, +1);
            cut_update_refs(ci, -1);

            if (ec < best_ec) {
                best_ec = ec; best_idx = ci; has_valid = true;
            }
        }

        if (has_valid && best_idx != old_idx) {
            selected_match_[n] = best_idx;
            if (ref_[n] > 0) cut_update_refs(best_idx, +1);
            auto const& best_match = cut_manager_.get_cut_matches()[best_idx];
            if (best_match.outputs.size() > 1) {
                for (xag_node o : best_match.outputs) {
                    if (o != n) {
                        selected_match_[o] = best_idx;
                        covered_.insert(o);
                        for (auto& u : best_match.children)
                            ref_[u] += 1;
                    }
                }
            }
        } else if (old_idx != 0 && ref_[n] > 0) {
            cut_update_refs(old_idx, +1);
        }
    }
}

uint32_t DaoMapOptimizer::cut_update_refs(uint32_t ci, int delta) {
    if (ci == 0 || ci >= cut_manager_.get_cut_matches().size()) return 0;
    auto& match = cut_manager_.get_cut_matches()[ci];
    auto sat_add = [](uint32_t a, uint32_t b) -> uint32_t {
        if (a == UINT32_MAX || b == UINT32_MAX) return UINT32_MAX;
        return (UINT32_MAX - a < b) ? UINT32_MAX : a + b;
    };
    uint32_t count = match.c.t_count;
    for (auto& child : match.children) {
        if (xag_.is_constant(child) || xag_.is_pi(child)) continue;
        bool should_recurse = (delta > 0) ? (ref_[child]++ == 0) : (--ref_[child] == 0);
        if (should_recurse && selected(child))
            count = sat_add(count, cut_update_refs(selected_match_[child], delta));
    }
    return count;
}

std::unordered_map<xag_node, uint32_t> DaoMapOptimizer::compute_mapped_refs() const {
    std::unordered_map<xag_node, uint32_t> refs;
    std::queue<xag_node> q;
    std::unordered_set<xag_node> visited;

    xag_.foreach_po([&](auto po) {
        xag_node pn = xag_.get_node(po);
        refs[pn]++; q.push(pn);
    });
    while (!q.empty()) {
        xag_node n = q.front(); q.pop();
        if (visited.count(n)) continue;
        visited.insert(n);
        if (xag_.is_constant(n) || xag_.is_pi(n)) continue;
        CutMatch const* match = selected(n);
        if (match) {
            if (match->outputs.size() > 1)
                for (xag_node o : match->outputs) visited.insert(o);
            for (auto const& child : match->children) {
                refs[child]++; q.push(child);
            }
        } else {
            xag_.foreach_fanin(n, [&](auto fi) {
                xag_node fn = xag_.get_node(fi);
                refs[fn]++; q.push(fn);
            });
        }
    }
    return refs;
}

double DaoMapOptimizer::compute_af_score_for_cut(CutMatch const& match, xag_node node) const {
    double ac = static_cast<double>(match.c.t_count);
    double ffs = 0.0;
    for (auto const& child : match.children) {
        if (xag_.is_constant(child) || xag_.is_pi(child)) continue;
        double af_child = (AF_[child] > 0.0) ? AF_[child] : params_.area_flow_init;
        ffs += af_child;
    }
    uint32_t ru = (ref_[node] > 0) ? ref_[node] : 1;
    return (ac + ffs) / compute_refs_eff(ru, need_[node], params_.alpha);
}

}  // namespace exact_t
