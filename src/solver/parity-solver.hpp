/*! \file parity-solver.hpp
 *  \brief A* search solver for parity network synthesis.
 */

#pragma once

#include "core/circuit.hpp"
#include "solver/graph.hpp"
#include "solver/parity-matrix.hpp"
#include "solver/solver.hpp"

#include <chrono>
#include <cstdint>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace exact_t {

/*! \brief Configuration for ParitySolver. */
struct ParitySolverParam {
    bool verbose{false};
    uint32_t n_branches{10};
    uint32_t max_iterations{1000000};
    std::string cache_path;
};

/*! \brief Statistics for ParitySolver performance. */
struct ParitySolverStats {
    uint64_t total_calls{0};
    uint64_t cache_lookups{0};
    uint64_t cache_hits{0};
    uint64_t cache_misses{0};
    double total_time_ms{0.0};
    double astar_time_ms{0.0};
    double gray_synth_time_ms{0.0};

    double cache_hit_rate() const {
        return cache_lookups > 0
            ? static_cast<double>(cache_hits) / cache_lookups : 0.0;
    }
};

/*! \brief A* search state for parity synthesis. */
struct PSearchState {
    ParityMatrix pmat;
    uint32_t cost{0};
    uint32_t heuristic{0};
    uint32_t total_cost() const { return cost + heuristic; }
    bool operator>(PSearchState const& other) const {
        return total_cost() > other.total_cost();
    }
};

/*! \brief Parent pointer for A* path reconstruction. */
struct ParentInfo {
    ParityMatrix parent_canonical;
    bool is_phase_gate{false};
    uint32_t control{0};
    uint32_t target{0};
    uint32_t phase_row{0};
    uint32_t phase_col{0};
    uint32_t dist_to_goal{0};
    std::vector<uint32_t> row_perm;
    std::vector<uint32_t> col_perm;
};

/*! \brief A* search solver for parity network synthesis. */
class ParitySolver {
  public:
    explicit ParitySolver(uint32_t n,
                          ParitySolverParam const& params = ParitySolverParam());

    /*! \brief Synthesize circuit for Z8 phase polynomial. */
    Circuit synthesize(Z8Phase const& z8_phase) const;

    uint32_t get_n() const { return n_; }
    ParitySolverStats const& get_stats() const { return stats_; }
    void reset_stats() { stats_ = ParitySolverStats(); }

    Graph& get_graph() { return graph_; }
    Graph const& get_graph() const { return graph_; }

    void load_cache(std::string const& path);
    void save_cache(std::string const& path) const;

  private:
    /*! \brief Mutable A* search context, allocated per astar_synth call. */
    struct AStarContext {
        std::priority_queue<PSearchState, std::vector<PSearchState>,
                            std::greater<PSearchState>> pq;
        std::unordered_set<ParityMatrix, ParityMatrixHash> visited;
        std::unordered_map<ParityMatrix, ParentInfo, ParityMatrixHash> parent_map;
    };

    Circuit gray_synth(ParityMatrix const& pmat) const;
    Circuit astar_synth(ParityMatrix const& pmat) const;

    bool process_phases(ParityMatrix& pmat, ParityMatrix& canon,
                        AStarContext& ctx) const;
    Circuit build_path(ParityMatrix const& goal_canon,
                       ParityMatrix const& initial_pmat,
                       AStarContext const& ctx) const;
    std::vector<std::tuple<uint32_t, uint32_t, int>>
    gen_cnot_candidates(ParityMatrix const& pmat, uint32_t n) const;
    void expand_cnot(ParityMatrix const& current_pmat,
                     ParityMatrix const& current_canonical,
                     uint32_t cost, uint32_t ctrl, uint32_t tgt,
                     AStarContext& ctx) const;

    uint32_t n_;
    ParitySolverParam params_;
    mutable Graph graph_;
    mutable ParitySolverStats stats_;
};

}  // namespace exact_t
