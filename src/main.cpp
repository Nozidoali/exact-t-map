#include "io/circuit-io.hpp"
#include "io/io.hpp"
#include "io/stats.hpp"
#include "mapper/mapper.hpp"

#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

struct Options {
    std::string input_file;
    std::string output_file;
    std::string format = "qasm";
    std::string cache_dir;
    uint32_t max_solver_n = 6;
    uint32_t cut_size = 5;
    uint32_t max_cuts_per_node = 20;
    uint32_t multioutput_cut_limit = 10;
    uint32_t area_flow_rounds = 6;
    uint32_t exact_rounds = 6;
    double alpha = 0.5;
    double beta = 0.1;
    double area_flow_init = 1.0;
    bool clean_ancilla = false;
    bool toffoli = false;
    bool multi_output = false;
    bool verbose = false;
    std::string stats_out;
    std::string trace_out;
    std::string dump_extraction;
};

void print_usage(char const* prog) {
    std::cerr << "Usage: " << prog << " [options] <input-file>\n"
              << "Options:\n"
              << "  -o <file>    Output file (default: stdout)\n"
              << "  -f qasm|qc   Output format (default: qasm)\n"
              << "  -n <int>     Max solver variables (default: 6)\n"
              << "  --cut-size <int>       Cut size (default: 5)\n"
              << "  --max-cuts-per-node <int>  Per-(node, size) cut cap (default: 20)\n"
              << "  --multioutput-cut-limit <int>  Per-(fanin, size) merge cap (default: 10)\n"
              << "  --area-flow-rounds <int>  Area flow rounds (default: 6)\n"
              << "  --exact-rounds <int>   Exact rounds (default: 6)\n"
              << "  --alpha <float>        Alpha (default: 0.5)\n"
              << "  --beta <float>         Beta (default: 0.1)\n"
              << "  --af-init <float>      Area flow init (default: 1.0)\n"
              << "  --clean-ancilla  Use clean ancilla don't-cares\n"
              << "  --toffoli    Use toffoli mapping\n"
              << "  --multi-output   Enable multi-output cut merging\n"
              << "  --verbose    Print statistics\n"
              << "  --stats-out <file>  Write RunStats JSON after mapping\n"
              << "  --trace-out <file>  Write per-round per-node JSONL trace\n"
              << "  --dump-extraction <file>  Write ILP extraction dump JSON after mapping\n";
}

bool parse_args(int argc, char* argv[], Options& opts) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            opts.output_file = argv[++i];
        } else if (std::strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            opts.format = argv[++i];
        } else if (std::strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            opts.max_solver_n = static_cast<uint32_t>(std::atoi(argv[++i]));
        } else if (std::strcmp(argv[i], "--cache-dir") == 0 && i + 1 < argc) {
            opts.cache_dir = argv[++i];
        } else if (std::strcmp(argv[i], "--cut-size") == 0 && i + 1 < argc) {
            opts.cut_size = static_cast<uint32_t>(std::atoi(argv[++i]));
        } else if (std::strcmp(argv[i], "--max-cuts-per-node") == 0 && i + 1 < argc) {
            opts.max_cuts_per_node = static_cast<uint32_t>(std::atoi(argv[++i]));
        } else if (std::strcmp(argv[i], "--multioutput-cut-limit") == 0 && i + 1 < argc) {
            opts.multioutput_cut_limit = static_cast<uint32_t>(std::atoi(argv[++i]));
        } else if (std::strcmp(argv[i], "--area-flow-rounds") == 0 && i + 1 < argc) {
            opts.area_flow_rounds = static_cast<uint32_t>(std::atoi(argv[++i]));
        } else if (std::strcmp(argv[i], "--exact-rounds") == 0 && i + 1 < argc) {
            opts.exact_rounds = static_cast<uint32_t>(std::atoi(argv[++i]));
        } else if (std::strcmp(argv[i], "--alpha") == 0 && i + 1 < argc) {
            opts.alpha = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "--beta") == 0 && i + 1 < argc) {
            opts.beta = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "--af-init") == 0 && i + 1 < argc) {
            opts.area_flow_init = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "--clean-ancilla") == 0) {
            opts.clean_ancilla = true;
        } else if (std::strcmp(argv[i], "--toffoli") == 0) {
            opts.toffoli = true;
        } else if (std::strcmp(argv[i], "--multi-output") == 0) {
            opts.multi_output = true;
        } else if (std::strcmp(argv[i], "--verbose") == 0) {
            opts.verbose = true;
        } else if (std::strcmp(argv[i], "--stats-out") == 0 && i + 1 < argc) {
            opts.stats_out = argv[++i];
        } else if (std::strcmp(argv[i], "--trace-out") == 0 && i + 1 < argc) {
            opts.trace_out = argv[++i];
        } else if (std::strcmp(argv[i], "--dump-extraction") == 0 && i + 1 < argc) {
            opts.dump_extraction = argv[++i];
        } else if (argv[i][0] == '-') {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            return false;
        } else {
            opts.input_file = argv[i];
        }
    }
    return !opts.input_file.empty();
}

}  // namespace

int main(int argc, char* argv[]) {
    Options opts;
    if (!parse_args(argc, argv, opts)) {
        print_usage(argv[0]);
        return 1;
    }

    auto xag = exact_t::read_network(opts.input_file);

    exact_t::MapperParams params;
    params.max_solver_n = opts.max_solver_n;
    params.cut_size = opts.cut_size;
    params.max_cuts_per_node = opts.max_cuts_per_node;
    params.multioutput_cut_limit = opts.multioutput_cut_limit;
    params.area_flow_rounds = opts.area_flow_rounds;
    params.exact_rounds = opts.exact_rounds;
    params.alpha = opts.alpha;
    params.beta = opts.beta;
    params.area_flow_init = opts.area_flow_init;
    params.use_clean_ancilla = opts.clean_ancilla;
    params.toffoli_mapping = opts.toffoli;
    params.multi_output = opts.multi_output;
    params.verbose = opts.verbose;
    params.stats_out_path = opts.stats_out;
    params.trace_out_path = opts.trace_out;
    params.dump_extraction_path = opts.dump_extraction;
    params.benchmark_name = std::filesystem::path(opts.input_file).filename().string();

    std::string cache_dir = opts.cache_dir;
    if (cache_dir.empty()) {
        auto candidate = std::filesystem::path(argv[0]).parent_path().parent_path() / "cache";
        if (std::filesystem::exists(candidate))
            cache_dir = candidate.string();
    }
    if (!cache_dir.empty()) {
        std::string solver_csv = cache_dir + "/solver/solver_cache_n" +
                                  std::to_string(opts.max_solver_n) + "_noancilla.csv";
        if (std::filesystem::exists(solver_csv))
            params.solver_cache_path = solver_csv;
        std::string parity_base = cache_dir + "/parity/parity_cache_n" +
                                   std::to_string(opts.max_solver_n) + "_noancilla_n" +
                                   std::to_string(opts.max_solver_n);
        if (std::filesystem::exists(parity_base + ".txt") ||
            std::filesystem::exists(parity_base + "_n" + std::to_string(opts.max_solver_n) + ".txt"))
            params.parity_cache_path = parity_base;
    }

    exact_t::Mapper mapper(xag, params);
    exact_t::Circuit circuit = mapper.map();

    if (opts.verbose)
        std::cerr << exact_t::format_run_stats_summary(mapper.get_run_stats());

    if (!opts.stats_out.empty())
        exact_t::write_run_stats_json(mapper.get_run_stats(), opts.stats_out);

    std::string output;
    if (opts.format == "qc")
        output = exact_t::to_qc(circuit);
    else
        output = exact_t::to_qasm(circuit);

    if (opts.output_file.empty()) {
        std::cout << output;
    } else {
        if (opts.format == "qc")
            exact_t::write_qc(circuit, opts.output_file);
        else
            exact_t::write_qasm(circuit, opts.output_file);
    }

    return 0;
}
