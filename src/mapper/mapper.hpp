/*! \file mapper.hpp
 *  \brief XAG to quantum circuit mapper using DAOMap optimization.
 */

#pragma once

#include <cstdint>

#include "io/stats.hpp"
#include "mapper/daomap.hpp"
#include "mapper/extractor.hpp"
#include "solver/parity-solver.hpp"

#include <chrono>
#include <string>

namespace exact_t {

/*! \brief Mapper configuration. */
struct MapperParams {
    uint32_t cut_size{5};
    uint32_t max_solver_n{6};
    uint32_t area_flow_rounds{6};
    uint32_t exact_rounds{6};
    double alpha{0.5};
    double beta{0.1};
    double area_flow_init{1.0};
    uint32_t t_count_per_and{7};
    uint32_t max_cuts_per_node{20};
    uint32_t multioutput_cut_limit{10};
    bool use_clean_ancilla{false};
    bool multi_output{false};
    bool toffoli_mapping{false};
    bool trivial_mapping{false};
    bool verbose{false};
    AllocationStrategy allocation{AllocationStrategy::LTFI};
    std::string solver_cache_path;
    std::string parity_cache_path;
    std::string stats_out_path;
    std::string trace_out_path;
    std::string dump_extraction_path;
    std::string benchmark_name;
};

/*! \brief XAG to quantum circuit mapper (thin orchestrator). */
class Mapper {
  public:
    Mapper(xag_network const& xag, MapperParams const& params = {});

    /*! \brief Execute mapping and return Clifford+T circuit. */
    Circuit map();

    /*! \brief Get mapping statistics. */
    MapperStats const& get_stats() const { return stats_; }

    /*! \brief Get aggregated run statistics (populated at end of map()). */
    RunStats const& get_run_stats() const { return run_stats_; }

    /*! \brief Access the cut manager (read-only). */
    CutManager const& cut_manager() const { return cut_manager_; }

    /*! \brief Access the XAG network (read-only). */
    xag_network const& xag() const { return xag_; }

  private:
    xag_network const& xag_;
    topo_view topo_;
    CutManager cut_manager_;
    MapperParams params_;
    CliffordTSolver solver_;
    ParitySolver parity_solver_;
    MapperStats stats_;
    RunStats run_stats_;

    void populate_run_stats(DaoMapOptimizer const& optimizer);
};

}  // namespace exact_t
