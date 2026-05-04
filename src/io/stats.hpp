/*! \file stats.hpp
 *  \brief Aggregated RunStats and its JSON / text emitters.
 */

#pragma once

#include "mapper/cut.hpp"
#include "mapper/mapper-stats.hpp"
#include "solver/parity-solver.hpp"
#include "solver/affine.hpp"

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>

namespace exact_t {

/*! \brief CliffordTSolver cache counters. */
struct SolverCacheStats {
    uint64_t size{0};
    uint64_t hits{0};
    uint64_t misses{0};
    double hit_rate{0.0};
};

/*! \brief Distribution data derived from final cut selections. */
struct DerivedStats {
    uint32_t self_cut_nodes{0};
    uint32_t mapped_and_nodes{0};
    double self_cut_ratio{0.0};
    std::map<uint32_t, uint32_t> cut_size_histogram;
    std::map<uint32_t, uint32_t> num_outputs_histogram;
    std::map<uint32_t, uint32_t> degree_histogram;
};

/*! \brief Subset of MapperParams serialized in RunStats for reproducibility. */
struct RunParams {
    uint32_t cut_size{0};
    uint32_t max_solver_n{0};
    uint32_t area_flow_rounds{0};
    uint32_t exact_rounds{0};
    double alpha{0.0};
    double beta{0.0};
    bool multi_output{false};
    bool clean_ancilla{false};
    bool toffoli{false};
};

/*! \brief Aggregated per-run statistics, serialized via write_run_stats_json. */
struct RunStats {
    uint32_t schema_version{1};
    std::string benchmark;
    RunParams params;
    MapperStats mapper;
    CutEnumStats cut_enum;
    SolverCacheStats solver_cache;
    ParitySolverStats parity_solver;
    /*! \brief Reserved; Affine is currently orphan code — fields are
     *  always zero until the synthesis path is wired to Mapper. */
    affine::AffineStats affine;
    DerivedStats derived;
};

/*! \brief Write RunStats to a JSON file. Aborts on open failure. */
void write_run_stats_json(RunStats const& stats, std::string const& path);

/*! \brief Human-readable multi-line summary (for stderr under --verbose). */
std::string format_run_stats_summary(RunStats const& stats);

}  // namespace exact_t
