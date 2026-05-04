/*! \file daomap.hpp
 *  \brief DAOMap cut selection optimizer for XAG networks.
 */

#pragma once

#include <cstdint>

#include "mapper/evaluator.hpp"
#include "mapper/cut.hpp"
#include "solver/solver.hpp"
#include "io/trace.hpp"

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace exact_t {

/*! \brief DAOMap cut selection optimizer.
 *
 *  Performs area-flow optimization and exact mapping rounds
 *  to select the best cut for each node in the XAG.
 */
class DaoMapOptimizer {
  public:
    struct Params {
        uint32_t area_flow_rounds{6};
        uint32_t exact_rounds{6};
        double alpha{0.5};
        double beta{0.1};
        double area_flow_init{1.0};
        uint32_t max_solver_n{6};
        uint32_t t_count_per_and{7};
        bool use_clean_ancilla{false};
        bool toffoli_mapping{false};
        bool verbose{false};
        std::string trace_out_path;
    };

    DaoMapOptimizer(xag_network const& xag, topo_view& topo,
                     CutManager& cut_manager, CliffordTSolver& solver,
                     Params const& params);

    /*! \brief Run cut selection optimization. */
    void optimize();

    /*! \brief Get selected cut for a node. */
    CutMatch const* selected(xag_node n) const;

    /*! \brief Compute total mapped cost across all selected cuts. */
    Cost compute_mapped_total_cost() const;

    /*! \brief Emit a JSONL phase_marker line. No-op when no sink. */
    void emit_phase_marker(std::string const& phase, std::string const& event,
                           double time_ms);

  private:
    void reset_state();
    /*! \brief Emit per-node events and a round_summary for a mapping round. */
    void emit_round_events(std::string const& phase, uint32_t round);

    void initialize_area_flow();
    void compute_area_flow();
    void initialize_po_references();
    void update_references_from_exact_mapping(uint32_t iteration);
    std::vector<xag_node> compute_reverse_topological_order() const;
    void evaluate_and_select_cuts(uint32_t iteration = 0);
    void evaluate_and_select_cuts_exact();
    uint32_t cut_update_refs(uint32_t cut_idx, int delta);
    void smooth_est_refs(double coef,
                          std::unordered_map<xag_node, uint32_t> const& exact_refs);
    std::unordered_map<xag_node, uint32_t> compute_mapped_refs() const;
    double compute_af_score_for_cut(CutMatch const& match, xag_node node) const;

    xag_network const& xag_;
    topo_view& topo_;
    CutManager& cut_manager_;
    CliffordTSolver& solver_;
    Params params_;
    CostEvaluator cost_evaluator_;

    node_map<uint32_t> selected_match_;
    node_map<double> AF_;
    node_map<uint32_t> ref_;
    node_map<uint32_t> est_ref_;
    node_map<uint32_t> occ_;
    node_map<uint32_t> need_;
    node_map<uint32_t> prev_selected_;
    bool have_prev_selected_{false};
    std::vector<xag_node> reverse_topo_order_;
    std::unordered_set<xag_node> covered_;
    std::unique_ptr<TraceSink> sink_;
};

}  // namespace exact_t
