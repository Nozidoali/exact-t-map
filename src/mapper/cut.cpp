#include "mapper/cut.hpp"

#include <algorithm>
#include <chrono>
#include <map>
#include <set>

namespace exact_t {

namespace {

template <typename ScoreFn>
std::vector<uint32_t> select_top_n(std::vector<uint32_t> const& candidates,
                                    uint32_t limit, ScoreFn score_fn) {
    std::vector<uint32_t> result = candidates;
    std::sort(result.begin(), result.end(),
              [&](uint32_t a, uint32_t b) { return score_fn(a) < score_fn(b); });
    if (result.size() > limit) result.resize(limit);
    return result;
}

template <typename ScoreFn>
std::vector<uint32_t> select_top_per_size(std::vector<uint32_t> const& source,
                                            std::vector<CutMatch> const& cut_matches,
                                            uint32_t per_size_limit,
                                            ScoreFn score_fn) {
    std::map<uint32_t, std::vector<uint32_t>> by_size;
    for (uint32_t idx : source) {
        uint32_t s = static_cast<uint32_t>(cut_matches[idx].children.size());
        by_size[s].push_back(idx);
    }
    std::vector<uint32_t> result;
    for (auto& [size, bucket] : by_size) {
        if (bucket.size() > per_size_limit) {
            std::sort(bucket.begin(), bucket.end(),
                      [&](uint32_t a, uint32_t b) {
                          return score_fn(a) < score_fn(b);
                      });
            bucket.resize(per_size_limit);
        }
        for (uint32_t idx : bucket) result.push_back(idx);
    }
    return result;
}

}  // namespace

CutMatch::CutMatch(std::vector<xag_node> const& children, xag_node output_node,
                     ANF2<xag_node> const& anf)
    : c({0, 0, 0, 0, 0}), children(children), outputs({output_node}), anfs({anf}) {
    std::sort(this->children.begin(), this->children.end());
}

CutMatch::CutMatch(std::vector<xag_node> const& children,
                     std::vector<xag_node> const& outputs,
                     std::vector<ANF2<xag_node>> const& anfs)
    : c({0, 0, 0, 0, 0}), children(children), outputs(outputs), anfs(anfs) {
    std::sort(this->children.begin(), this->children.end());
}

std::string CutMatch::to_string() const {
    std::string result =
        "CutMatch(degree=" + std::to_string(degree()) +
        ", cost={T:" + std::to_string(c.t_count) +
        ", CX:" + std::to_string(c.cnot_count) +
        ", Q:" + std::to_string(c.q_count) + "}, children=[";
    for (size_t i = 0; i < children.size(); ++i) {
        if (i > 0) result += ", ";
        result += std::to_string(children[i]);
    }
    result += "], outputs=[";
    for (size_t i = 0; i < outputs.size(); ++i) {
        if (i > 0) result += ", ";
        result += std::to_string(outputs[i]);
    }
    result += "])";
    return result;
}

CutManager::CutManager(xag_network const& xag, topo_view const& topo,
                        CutEnumParam const& params)
    : xag_(xag), topo_(topo), params_(params), candidate_matches_(xag) {}

void CutManager::cut_enum_bilinear() {
    auto start = std::chrono::high_resolution_clock::now();
    stats_.reset();
    cut_matches_.clear();
    topo_.foreach_node([&](auto n) { candidate_matches_[n].clear(); });

    topo_.foreach_node([&](auto n) {
        add_match(CutMatch({n}, n, ANF2<xag_node>(n)), n);

        if (!xag_.is_and(n) && !xag_.is_xor(n)) return;
        std::vector<xag_signal> fanins;
        xag_.foreach_fanin(n, [&](auto fi) { fanins.push_back(fi); });
        if (fanins.size() != 2) return;
        bool is_and = xag_.is_and(n);
        xag_node const f0 = xag_.get_node(fanins[0]);
        xag_node const f1 = xag_.get_node(fanins[1]);

        auto find_self_cut = [&](xag_node fn) -> uint32_t {
            for (uint32_t idx : candidate_matches_[fn]) {
                CutMatch const& m = cut_matches_[idx];
                if (m.is_self_cut() && !m.outputs.empty() && m.outputs[0] == fn)
                    return idx;
            }
            return 0;
        };
        {
            uint32_t s0 = find_self_cut(f0);
            uint32_t s1 = find_self_cut(f1);
            if (s0 != 0 && s1 != 0) {
                auto merged = merge_children(cut_matches_[s0],
                                              cut_matches_[s1], params_.max_cut_size);
                if (!merged.empty()) {
                    merge_cuts_impl(n, s0, s1, merged,
                                    xag_.is_complemented(fanins[0]),
                                    xag_.is_complemented(fanins[1]), is_and);
                    stats_.trivial_cuts++;
                }
            }
        }

        if (params_.only_trivial_cuts) return;

        auto cuts_f0 = select_top_candidate_matches(candidate_matches_[f0]);
        auto cuts_f1 = select_top_candidate_matches(candidate_matches_[f1]);

        for (uint32_t cut0_idx : cuts_f0) {
            for (uint32_t cut1_idx : cuts_f1) {
                CutMatch const& match0 = cut_matches_[cut0_idx];
                CutMatch const& match1 = cut_matches_[cut1_idx];
                if (is_and && (match0.degree() > 1 || match1.degree() > 1))
                    continue;
                auto merged = merge_children(match0, match1, params_.max_cut_size);
                if (!merged.empty()) {
                    merge_cuts_impl(n, cut0_idx, cut1_idx, merged,
                                    xag_.is_complemented(fanins[0]),
                                    xag_.is_complemented(fanins[1]), is_and);
                    stats_.bilinear_cuts++;
                }
            }
        }
        prune_node_cuts_by_score(n);
    });

    topo_.foreach_node([&](auto n) {
        auto& matches = candidate_matches_[n];
        matches.erase(std::remove_if(matches.begin(), matches.end(),
                                     [&](uint32_t idx) {
                                         auto const& m = cut_matches_[idx];
                                         return m.degree() == 0 &&
                                                m.children.size() == 1 &&
                                                m.children[0] == n;
                                     }),
                      matches.end());
    });

    if (params_.enable_multi_output_cells) add_multi_output_cuts();

    topo_.foreach_node([&](auto n) {
        if (!xag_.is_constant(n) && !xag_.is_pi(n)) {
            auto& matches = candidate_matches_[n];
            stats_.cuts_by_output_count[static_cast<uint32_t>(matches.size())]++;
            matches = select_top_per_size(matches, cut_matches_,
                                           params_.max_cuts_per_node,
                                           [&](uint32_t idx) {
                                               return get_cut_candidate_priority(idx);
                                           });
        }
    });

    stats_.total_cuts = cut_matches_.size();
    stats_.avg_cuts_per_node =
        static_cast<double>(stats_.total_cuts) / std::max(1u, xag_.num_gates());

    auto end = std::chrono::high_resolution_clock::now();
    stats_.time_build_cuts_ms =
        std::chrono::duration<double, std::milli>(end - start).count();
}

uint32_t CutManager::add_match(CutMatch const& match, xag_node node) {
    for (uint32_t idx : candidate_matches_[node])
        if (cut_matches_[idx] == match) return idx;
    uint32_t index = static_cast<uint32_t>(cut_matches_.size());
    cut_matches_.push_back(match);
    candidate_matches_[node].push_back(index);
    return index;
}

uint32_t CutManager::add_match_and_prune(CutMatch const& match, xag_node node) {
    uint32_t idx = add_match(match, node);
    prune_node_cuts_by_score(node);
    return idx;
}

uint32_t CutManager::add_multi_output_match(CutMatch const& match) {
    for (auto out : match.outputs) {
        for (uint32_t idx : candidate_matches_[out])
            if (cut_matches_[idx] == match) return idx;
    }
    uint32_t index = static_cast<uint32_t>(cut_matches_.size());
    cut_matches_.push_back(match);
    for (auto out : match.outputs)
        candidate_matches_[out].push_back(index);
    return index;
}

void CutManager::add_multi_output_cuts() {
    struct SiblingEntry {
        xag_node and_node{0};
        xag_node xor_node{0};
        bool and_compl0{false};
        bool and_compl1{false};
        bool xor_compl0{false};
        bool xor_compl1{false};
    };

    std::map<std::pair<xag_node, xag_node>, SiblingEntry> sibling_map;

    topo_.foreach_node([&](auto n) {
        if (xag_.is_constant(n) || xag_.is_pi(n)) return;
        if (!xag_.is_and(n) && !xag_.is_xor(n)) return;
        std::vector<xag_signal> fanins;
        xag_.foreach_fanin(n, [&](auto fi) { fanins.push_back(fi); });
        if (fanins.size() != 2) return;

        xag_node f0 = xag_.get_node(fanins[0]);
        xag_node f1 = xag_.get_node(fanins[1]);
        xag_node lo = std::min(f0, f1);
        xag_node hi = std::max(f0, f1);
        auto key = std::make_pair(lo, hi);

        bool c0 = xag_.is_complemented(fanins[0]);
        bool c1 = xag_.is_complemented(fanins[1]);
        if (f0 > f1) { std::swap(c0, c1); }

        auto& entry = sibling_map[key];
        if (xag_.is_and(n)) {
            entry.and_node = n;
            entry.and_compl0 = c0;
            entry.and_compl1 = c1;
        } else {
            entry.xor_node = n;
            entry.xor_compl0 = c0;
            entry.xor_compl1 = c1;
        }
    });

    for (auto const& [key, entry] : sibling_map) {
        if (entry.and_node == 0 || entry.xor_node == 0) continue;

        auto cuts_f0 = select_top_candidate_matches(candidate_matches_[key.first]);
        auto cuts_f1 = select_top_candidate_matches(candidate_matches_[key.second]);

        for (uint32_t c0 : cuts_f0) {
            for (uint32_t c1 : cuts_f1) {
                auto const& m0 = cut_matches_[c0];
                auto const& m1 = cut_matches_[c1];
                if (m0.degree() > 1 || m1.degree() > 1) continue;

                auto merged = merge_children(m0, m1);
                if (merged.size() > params_.max_cut_size) continue;

                ANF2<xag_node> anf0 = m0.anfs[0];
                ANF2<xag_node> anf1 = m1.anfs[0];

                ANF2<xag_node> and_anf = ANF2<xag_node>::merge(
                    anf0, anf1, entry.and_compl0, entry.and_compl1, true);
                ANF2<xag_node> xor_anf = ANF2<xag_node>::merge(
                    anf0, anf1, entry.xor_compl0, entry.xor_compl1, false);

                std::set<xag_node> anf_vars;
                for (auto const& t : and_anf.degree1) anf_vars.insert(t);
                for (auto const& [i, j] : and_anf.degree2) {
                    anf_vars.insert(i);
                    anf_vars.insert(j);
                }
                for (auto const& t : xor_anf.degree1) anf_vars.insert(t);
                for (auto const& [i, j] : xor_anf.degree2) {
                    anf_vars.insert(i);
                    anf_vars.insert(j);
                }

                std::vector<xag_node> actual(anf_vars.begin(), anf_vars.end());
                if (actual.size() > params_.max_cut_size) continue;

                CutMatch match(actual,
                               {entry.and_node, entry.xor_node},
                               {and_anf, xor_anf});
                add_multi_output_match(match);
                stats_.multi_output_cuts++;
            }
        }
    }
}

std::vector<xag_node> CutManager::merge_children(CutMatch const& match0,
                                                   CutMatch const& match1,
                                                   uint32_t max_children) const {
    std::set<xag_node> unique;
    for (auto n : match0.children) unique.insert(n);
    for (auto n : match1.children) unique.insert(n);
    std::vector<xag_node> merged;
    for (auto n : unique) {
        if (max_children > 0 && merged.size() >= max_children) break;
        merged.push_back(n);
    }
    return merged;
}

void CutManager::merge_cuts_impl(xag_node node, uint32_t cut0_idx, uint32_t cut1_idx,
                                  std::vector<xag_node> const& merged_children,
                                  bool compl0, bool compl1, bool is_and) {
    if (merged_children.size() > params_.max_cut_size) return;
    ANF2<xag_node> anf0 = cut_matches_[cut0_idx].anfs[0];
    ANF2<xag_node> anf1 = cut_matches_[cut1_idx].anfs[0];
    ANF2<xag_node> result_anf = ANF2<xag_node>::merge(anf0, anf1, compl0, compl1, is_and);

    std::set<xag_node> anf_nodes;
    for (auto const& t : result_anf.degree1) anf_nodes.insert(t);
    for (auto const& [i, j] : result_anf.degree2) {
        anf_nodes.insert(i);
        anf_nodes.insert(j);
    }

    std::vector<xag_node> actual(anf_nodes.begin(), anf_nodes.end());
    if (actual.size() > params_.max_cut_size) return;

    CutMatch match(actual, node, result_anf);
    add_match_and_prune(match, node);
    stats_.merged_cuts++;
}

void CutManager::prune_node_cuts_by_score(xag_node node) {
    auto& matches = candidate_matches_[node];
    matches = select_top_per_size(matches, cut_matches_,
                                   params_.max_cuts_per_node_during_merge,
                                   [&](uint32_t idx) {
                                       auto const& m = cut_matches_[idx];
                                       return m.score > 0.0
                                           ? m.score
                                           : static_cast<double>(m.c.t_count);
                                   });
}

std::vector<uint32_t> CutManager::select_top_candidate_matches(
    std::vector<uint32_t> const& source) const {
    return select_top_per_size(source, cut_matches_, params_.multioutput_cut_limit,
                               [&](uint32_t idx) { return get_cut_candidate_priority(idx); });
}

double CutManager::get_cut_candidate_priority(uint32_t cut_idx) const {
    auto const& m = cut_matches_[cut_idx];
    return static_cast<double>(m.c.t_count) * 1e9 +
           static_cast<double>(m.c.cnot_count) * 1e3 +
           static_cast<double>(m.c.q_count);
}

}  // namespace exact_t
