#include "io/stats.hpp"

#include <cassert>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace exact_t {

namespace {

std::string quote(std::string const& s) {
    std::string out = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') out += '\\';
        out += c;
    }
    out += '"';
    return out;
}

std::string fmt_double(double x) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << x;
    return oss.str();
}

std::string emit_histogram(std::map<uint32_t, uint32_t> const& h) {
    std::ostringstream oss;
    oss << '{';
    bool first = true;
    for (auto const& [k, v] : h) {
        if (!first) oss << ',';
        first = false;
        oss << '"' << k << "\":" << v;
    }
    oss << '}';
    return oss.str();
}

}  // namespace

void write_run_stats_json(RunStats const& s, std::string const& path) {
    std::ofstream out(path);
    assert(out.is_open());
    out << std::fixed << std::setprecision(6);

    out << "{\n";
    out << "  \"schema_version\":" << s.schema_version << ",\n";
    out << "  \"benchmark\":" << quote(s.benchmark) << ",\n";

    out << "  \"params\":{";
    out << "\"cut_size\":" << s.params.cut_size
        << ",\"max_solver_n\":" << s.params.max_solver_n
        << ",\"area_flow_rounds\":" << s.params.area_flow_rounds
        << ",\"exact_rounds\":" << s.params.exact_rounds
        << ",\"alpha\":" << fmt_double(s.params.alpha)
        << ",\"beta\":" << fmt_double(s.params.beta)
        << ",\"multi_output\":" << (s.params.multi_output ? "true" : "false")
        << ",\"clean_ancilla\":" << (s.params.clean_ancilla ? "true" : "false")
        << ",\"toffoli\":" << (s.params.toffoli ? "true" : "false")
        << "},\n";

    out << "  \"mapper\":{"
        << "\"baseline_and_gates\":" << s.mapper.baseline_and_gates
        << ",\"baseline_t_count\":" << s.mapper.baseline_t_count
        << ",\"mapped_t_count\":" << s.mapper.mapped_t_count
        << ",\"mapped_cnot_count\":" << s.mapper.mapped_cnot_count
        << ",\"mapped_qubit_count\":" << s.mapper.mapped_qubit_count
        << ",\"runtime_ms\":" << fmt_double(s.mapper.runtime_ms)
        << "},\n";

    out << "  \"cut_enum\":{"
        << "\"total_cuts\":" << s.cut_enum.total_cuts
        << ",\"trivial_cuts\":" << s.cut_enum.trivial_cuts
        << ",\"bilinear_cuts\":" << s.cut_enum.bilinear_cuts
        << ",\"multi_output_cuts\":" << s.cut_enum.multi_output_cuts
        << ",\"merged_cuts\":" << s.cut_enum.merged_cuts
        << ",\"avg_cuts_per_node\":" << fmt_double(s.cut_enum.avg_cuts_per_node)
        << ",\"time_build_cuts_ms\":" << fmt_double(s.cut_enum.time_build_cuts_ms)
        << ",\"cuts_by_output_count\":" << emit_histogram(s.cut_enum.cuts_by_output_count)
        << "},\n";

    out << "  \"solver_cache\":{"
        << "\"size\":" << s.solver_cache.size
        << ",\"hits\":" << s.solver_cache.hits
        << ",\"misses\":" << s.solver_cache.misses
        << ",\"hit_rate\":" << fmt_double(s.solver_cache.hit_rate)
        << "},\n";

    out << "  \"parity_solver\":{"
        << "\"total_calls\":" << s.parity_solver.total_calls
        << ",\"cache_lookups\":" << s.parity_solver.cache_lookups
        << ",\"cache_hits\":" << s.parity_solver.cache_hits
        << ",\"cache_misses\":" << s.parity_solver.cache_misses
        << ",\"cache_hit_rate\":" << fmt_double(s.parity_solver.cache_hit_rate())
        << ",\"total_time_ms\":" << fmt_double(s.parity_solver.total_time_ms)
        << ",\"astar_time_ms\":" << fmt_double(s.parity_solver.astar_time_ms)
        << ",\"gray_synth_time_ms\":" << fmt_double(s.parity_solver.gray_synth_time_ms)
        << "},\n";

    out << "  \"affine\":{"
        << "\"initial_cnot_count\":" << s.affine.initial_cnot_count
        << ",\"final_cnot_count\":" << s.affine.final_cnot_count
        << ",\"blocks_processed\":" << s.affine.blocks_processed
        << ",\"states_explored\":" << s.affine.states_explored
        << ",\"timeouts\":" << s.affine.timeouts
        << ",\"total_time_ms\":" << fmt_double(s.affine.total_time_ms)
        << "},\n";

    out << "  \"derived\":{"
        << "\"self_cut_nodes\":" << s.derived.self_cut_nodes
        << ",\"mapped_and_nodes\":" << s.derived.mapped_and_nodes
        << ",\"self_cut_ratio\":" << fmt_double(s.derived.self_cut_ratio)
        << ",\"cut_size_histogram\":" << emit_histogram(s.derived.cut_size_histogram)
        << ",\"num_outputs_histogram\":" << emit_histogram(s.derived.num_outputs_histogram)
        << ",\"degree_histogram\":" << emit_histogram(s.derived.degree_histogram)
        << "}\n";

    out << "}\n";
}

namespace {

std::string format_histogram(std::map<uint32_t, uint32_t> const& h) {
    std::ostringstream oss;
    bool first = true;
    for (auto const& [k, v] : h) {
        if (!first) oss << "  ";
        first = false;
        oss << k << ":" << v;
    }
    return oss.str();
}

}  // namespace

std::string format_run_stats_summary(RunStats const& s) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);

    oss << "[params]    cut_size=" << s.params.cut_size
        << " n=" << s.params.max_solver_n
        << " af_rounds=" << s.params.area_flow_rounds
        << " exact_rounds=" << s.params.exact_rounds
        << " multi_output=" << (s.params.multi_output ? "on" : "off") << "\n";
    oss << "\n";

    oss << "[baseline]  AND=" << s.mapper.baseline_and_gates
        << "  T=" << s.mapper.baseline_t_count << "\n";

    int32_t dt = static_cast<int32_t>(s.mapper.mapped_t_count) -
                  static_cast<int32_t>(s.mapper.baseline_t_count);
    double dt_pct = s.mapper.baseline_t_count > 0
        ? 100.0 * dt / static_cast<double>(s.mapper.baseline_t_count) : 0.0;
    oss << "[mapped]    T=" << s.mapper.mapped_t_count
        << "  CNOT=" << s.mapper.mapped_cnot_count
        << "  qubits=" << s.mapper.mapped_qubit_count
        << "  runtime=" << s.mapper.runtime_ms << "ms"
        << "  (Δt=" << dt << ", " << std::setprecision(2) << dt_pct << "%)\n";
    oss << "\n";

    oss << "[cut_enum]  total=" << s.cut_enum.total_cuts
        << "  trivial=" << s.cut_enum.trivial_cuts
        << "  bilinear=" << s.cut_enum.bilinear_cuts
        << "  multi=" << s.cut_enum.multi_output_cuts
        << "  avg/node=" << s.cut_enum.avg_cuts_per_node
        << "  time=" << s.cut_enum.time_build_cuts_ms << "ms\n";

    oss << "[selected]  self_cut=" << s.derived.self_cut_nodes
        << " (" << std::setprecision(1) << (s.derived.self_cut_ratio * 100.0) << "%)"
        << "  mapped_and=" << s.derived.mapped_and_nodes
        << "/" << s.mapper.baseline_and_gates << "\n";
    oss << "            size hist:   " << format_histogram(s.derived.cut_size_histogram) << "\n";
    oss << "            degree hist: " << format_histogram(s.derived.degree_histogram) << "\n";
    oss << "\n";

    oss << std::setprecision(2);
    oss << "[solver]    cache hits=" << s.solver_cache.hits
        << " misses=" << s.solver_cache.misses
        << " hit_rate=" << (s.solver_cache.hit_rate * 100.0) << "%\n";

    oss << "[parity]    calls=" << s.parity_solver.total_calls
        << " hit_rate=" << (s.parity_solver.cache_hit_rate() * 100.0) << "%"
        << " time=" << s.parity_solver.total_time_ms << "ms"
        << " (astar=" << s.parity_solver.astar_time_ms << "ms, "
        << "gray=" << s.parity_solver.gray_synth_time_ms << "ms)\n";

    oss << "[affine] blocks=" << s.affine.blocks_processed
        << " states=" << s.affine.states_explored
        << " timeouts=" << s.affine.timeouts
        << " time=" << s.affine.total_time_ms << "ms\n";

    return oss.str();
}

}  // namespace exact_t
