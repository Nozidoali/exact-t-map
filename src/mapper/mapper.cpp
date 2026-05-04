#include "mapper/mapper.hpp"

#include <chrono>

#include "io/extraction-dump.hpp"

namespace exact_t {

Mapper::Mapper(xag_network const& xag, MapperParams const& params)
    : xag_(xag), topo_(xag_),
      cut_manager_(xag, topo_,
                   [&]() {
                       CutEnumParam cp;
                       cp.max_cut_size = params.cut_size;
                       cp.max_cuts_per_node = params.max_cuts_per_node;
                       cp.max_cuts_per_node_during_merge = params.max_cuts_per_node;
                       cp.multioutput_cut_limit = params.multioutput_cut_limit;
                       cp.t_count_per_and = params.t_count_per_and;
                       cp.verbose = params.verbose;
                       cp.only_trivial_cuts = params.trivial_mapping;
                       cp.enable_multi_output_cells = params.multi_output;
                       return cp;
                   }()),
      params_(params),
      solver_(params.max_solver_n),
      parity_solver_(params.max_solver_n,
                     [&]() {
                         ParitySolverParam pp;
                         pp.verbose = params.verbose;
                         pp.cache_path = params.parity_cache_path;
                         return pp;
                     }()) {}

Circuit Mapper::map() {
    auto t_start = std::chrono::high_resolution_clock::now();

    DaoMapOptimizer::Params dao_params;
    dao_params.area_flow_rounds = params_.area_flow_rounds;
    dao_params.exact_rounds = params_.exact_rounds;
    dao_params.alpha = params_.alpha;
    dao_params.beta = params_.beta;
    dao_params.area_flow_init = params_.area_flow_init;
    dao_params.max_solver_n = params_.max_solver_n;
    dao_params.t_count_per_and = params_.t_count_per_and;
    dao_params.use_clean_ancilla = params_.use_clean_ancilla;
    dao_params.toffoli_mapping = params_.toffoli_mapping;
    dao_params.verbose = params_.verbose;
    dao_params.trace_out_path = params_.trace_out_path;

    DaoMapOptimizer optimizer(xag_, topo_, cut_manager_, solver_, dao_params);
    optimizer.optimize();

    CircuitExtractor extractor(xag_, optimizer, solver_, params_.toffoli_mapping,
                               params_.use_clean_ancilla);
    auto t_extract = std::chrono::high_resolution_clock::now();
    Circuit result = extractor.extract(params_.allocation);
    double extract_ms = std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - t_extract).count();
    optimizer.emit_phase_marker("extract", "done", extract_ms);

    auto t_end = std::chrono::high_resolution_clock::now();
    stats_.runtime_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    Cost tc = optimizer.compute_mapped_total_cost();
    stats_.mapped_t_count = tc.t_count;
    stats_.mapped_cnot_count = tc.cnot_count;
    stats_.mapped_qubit_count = tc.q_count;

    uint32_t and_count = 0;
    xag_.foreach_node([&](auto n) { if (xag_.is_and(n)) and_count++; });
    stats_.baseline_and_gates = and_count;
    stats_.baseline_t_count = and_count * params_.t_count_per_and;

    populate_run_stats(optimizer);

    if (!params_.dump_extraction_path.empty())
        write_extraction_dump(*this, optimizer, params_.dump_extraction_path);

    return result;
}

void Mapper::populate_run_stats(DaoMapOptimizer const& optimizer) {
    run_stats_.benchmark = params_.benchmark_name;
    run_stats_.params.cut_size = params_.cut_size;
    run_stats_.params.max_solver_n = params_.max_solver_n;
    run_stats_.params.area_flow_rounds = params_.area_flow_rounds;
    run_stats_.params.exact_rounds = params_.exact_rounds;
    run_stats_.params.alpha = params_.alpha;
    run_stats_.params.beta = params_.beta;
    run_stats_.params.multi_output = params_.multi_output;
    run_stats_.params.clean_ancilla = params_.use_clean_ancilla;
    run_stats_.params.toffoli = params_.toffoli_mapping;

    run_stats_.mapper = stats_;
    run_stats_.cut_enum = cut_manager_.get_stats();
    run_stats_.parity_solver = parity_solver_.get_stats();

    run_stats_.solver_cache.size = solver_.cache_size();
    run_stats_.solver_cache.hits = solver_.cache_hits();
    run_stats_.solver_cache.misses = solver_.cache_misses();
    uint64_t lookups = run_stats_.solver_cache.hits + run_stats_.solver_cache.misses;
    run_stats_.solver_cache.hit_rate = lookups > 0
        ? static_cast<double>(run_stats_.solver_cache.hits) / lookups : 0.0;

    run_stats_.derived.self_cut_nodes = 0;
    uint32_t mapped_nodes = 0;
    xag_.foreach_node([&](auto n) {
        if (xag_.is_constant(n) || xag_.is_pi(n) || !xag_.is_and(n)) return;
        CutMatch const* m = optimizer.selected(n);
        if (m == nullptr) return;
        ++mapped_nodes;
        if (m->is_self_cut()) ++run_stats_.derived.self_cut_nodes;
        run_stats_.derived.cut_size_histogram[
            static_cast<uint32_t>(m->children.size())]++;
        run_stats_.derived.num_outputs_histogram[
            static_cast<uint32_t>(m->outputs.size())]++;
        run_stats_.derived.degree_histogram[m->degree()]++;
    });
    run_stats_.derived.mapped_and_nodes = mapped_nodes;
    run_stats_.derived.self_cut_ratio = mapped_nodes > 0
        ? static_cast<double>(run_stats_.derived.self_cut_nodes) / mapped_nodes
        : 0.0;
}

}  // namespace exact_t
