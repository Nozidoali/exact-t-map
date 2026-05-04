/*! \file cut.hpp
 *  \brief Cut enumeration for XAG to quantum circuit mapping.
 */

#pragma once

#include "core/anf.hpp"
#include "network/xag.hpp"

#include <algorithm>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace exact_t {

/*! \brief Compute effective reference count. */
inline double compute_refs_eff(uint32_t refs, uint32_t need, double alpha) {
    double r = static_cast<double>(std::max(1u, refs)) +
               alpha * static_cast<double>(need);
    return r < 1e-9 ? 1e-9 : r;
}

/*! \brief Resource cost for a quantum circuit or cut. */
struct Cost {
    uint32_t t_count{0};
    uint32_t cnot_count{0};
    uint32_t q_count{0};
    uint32_t s_count{0};
    uint32_t z_count{0};

    Cost operator+(Cost const& o) const {
        return {t_count + o.t_count, cnot_count + o.cnot_count,
                q_count + o.q_count, s_count + o.s_count, z_count + o.z_count};
    }

    bool operator<(Cost const& o) const {
        if (t_count != o.t_count) return t_count < o.t_count;
        if (cnot_count != o.cnot_count) return cnot_count < o.cnot_count;
        if (q_count != o.q_count) return q_count < o.q_count;
        if (s_count != o.s_count) return s_count < o.s_count;
        return z_count < o.z_count;
    }
};

/*! \brief Cut match representing a subgraph for mapping. */
struct CutMatch {
    Cost c;
    double score{0.0};
    std::vector<xag_node> children;
    std::vector<xag_node> outputs;
    std::vector<ANF2<xag_node>> anfs;

    CutMatch(std::vector<xag_node> const& children, xag_node output_node,
              ANF2<xag_node> const& anf);

    CutMatch(std::vector<xag_node> const& children,
             std::vector<xag_node> const& outputs,
             std::vector<ANF2<xag_node>> const& anfs);

    /*! \brief Get maximum degree across all ANFs. */
    uint32_t degree() const {
        uint32_t d = 0;
        for (auto const& anf : anfs) d = std::max(d, anf.get_degree());
        return d;
    }

    bool operator==(CutMatch const& other) const {
        return children == other.children && outputs == other.outputs;
    }

    /*! \brief Hash functor. */
    struct Hash {
        size_t operator()(CutMatch const& m) const {
            size_t h = 0x9e3779b97f4a7c15ULL;
            for (auto const& node : m.children)
                h ^= std::hash<uint32_t>{}(node) + 0x9e37 + (h << 6) + (h >> 2);
            for (auto const& node : m.outputs)
                h ^= std::hash<uint32_t>{}(node) + 0x85eb + (h << 6) + (h >> 2);
            return h;
        }
    };

    /*! \brief Check if cut is trivial (single node). */
    bool is_self_cut() const {
        return children.size() == 1 && outputs.size() == 1 &&
               children[0] == outputs[0];
    }

    std::string to_string() const;
};

/*! \brief Parameters for cut enumeration. */
struct CutEnumParam {
    uint32_t max_cut_size{5};
    uint32_t min_cut_size{2};
    /*! \brief Per-(node, cut size) cap for retained cuts after enumeration.
     *  Total per-node cap is max_cuts_per_node × #sizes. */
    uint32_t max_cuts_per_node{20};
    /*! \brief Per-(node, cut size) cap during incremental merge pruning.
     *  Total per-node cap during merge is max_cuts_per_node_during_merge × #sizes. */
    uint32_t max_cuts_per_node_during_merge{20};
    bool enable_multi_output_cells{false};
    uint32_t multioutput_cut_limit{10};
    uint32_t t_count_per_and{7};
    bool use_clean_ancilla{false};
    bool verbose{false};
    bool only_trivial_cuts{false};
};

/*! \brief Statistics from cut enumeration. */
struct CutEnumStats {
    uint64_t total_cuts{0};
    uint64_t trivial_cuts{0};
    uint64_t bilinear_cuts{0};
    uint64_t multi_output_cuts{0};
    uint64_t merged_cuts{0};
    double avg_cuts_per_node{0.0};
    double time_build_cuts_ms{0.0};
    std::map<uint32_t, uint32_t> cuts_by_output_count;

    void reset() {
        total_cuts = trivial_cuts = bilinear_cuts = 0;
        multi_output_cuts = merged_cuts = 0;
        avg_cuts_per_node = time_build_cuts_ms = 0.0;
        cuts_by_output_count.clear();
    }
};

/*! \brief Manager for cut enumeration on XAG networks. */
class CutManager {
  public:
    CutManager(xag_network const& xag, topo_view const& topo,
               CutEnumParam const& params);

    /*! \brief Enumerate bilinear cuts on the network. */
    void cut_enum_bilinear();

    /*! \brief Add multi-output cuts by merging. */
    void add_multi_output_cuts();

    void merge_cuts_impl(xag_node node, uint32_t cut0_idx, uint32_t cut1_idx,
                         std::vector<xag_node> const& merged_children,
                         bool compl0, bool compl1, bool is_and);

    void prune_node_cuts_by_score(xag_node node);

    uint32_t add_match(CutMatch const& match, xag_node node);
    uint32_t add_match_and_prune(CutMatch const& match, xag_node node);
    uint32_t add_multi_output_match(CutMatch const& match);

    std::vector<xag_node> merge_children(CutMatch const& match0,
                                          CutMatch const& match1,
                                          uint32_t max_children = 0) const;

    std::vector<CutMatch> const& get_cut_matches() const { return cut_matches_; }
    std::vector<CutMatch>& get_cut_matches() { return cut_matches_; }
    node_map<std::vector<uint32_t>> const& get_candidate_matches() const {
        return candidate_matches_;
    }
    CutEnumStats const& get_stats() const { return stats_; }
    CutEnumStats& get_stats() { return stats_; }

  private:
    xag_network const& xag_;
    topo_view const& topo_;
    CutEnumParam params_;
    CutEnumStats stats_;

    std::vector<CutMatch> cut_matches_;
    node_map<std::vector<uint32_t>> candidate_matches_;

    double get_cut_candidate_priority(uint32_t cut_idx) const;
    std::vector<uint32_t> select_top_candidate_matches(
        std::vector<uint32_t> const& source) const;
};

}  // namespace exact_t
