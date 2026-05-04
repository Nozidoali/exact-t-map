#include "io/extraction-dump.hpp"

#include "mapper/cut.hpp"
#include "mapper/daomap.hpp"
#include "mapper/mapper.hpp"
#include "network/xag.hpp"

#include <cassert>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace exact_t {

namespace {

void emit_int_array(std::ostream& out, std::vector<uint32_t> const& v) {
    out << '[';
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) out << ',';
        out << v[i];
    }
    out << ']';
}

}  // namespace

void write_extraction_dump(Mapper const& mapper, DaoMapOptimizer const& optimizer,
                           std::string const& path) {
    std::ofstream out(path);
    assert(out.is_open());

    auto const& xag = mapper.xag();
    auto const& cut_manager = mapper.cut_manager();
    auto const& cut_matches = cut_manager.get_cut_matches();
    auto const& candidates = cut_manager.get_candidate_matches();
    auto const& run_stats = mapper.get_run_stats();

    std::vector<uint32_t> po_nodes;
    xag.foreach_po([&](auto po) { po_nodes.push_back(xag.get_node(po)); });

    std::vector<uint32_t> and_nodes, xor_nodes;
    xag.foreach_node([&](auto n) {
        if (xag.is_constant(n) || xag.is_pi(n)) return;
        if (xag.is_and(n)) and_nodes.push_back(n);
        else xor_nodes.push_back(n);
    });

    out << "{\n";
    out << "  \"schema_version\":1,\n";
    out << "  \"benchmark\":\"" << run_stats.benchmark << "\",\n";

    out << "  \"po_nodes\":"; emit_int_array(out, po_nodes); out << ",\n";
    out << "  \"and_nodes\":"; emit_int_array(out, and_nodes); out << ",\n";
    out << "  \"xor_nodes\":"; emit_int_array(out, xor_nodes); out << ",\n";

    out << "  \"nodes\":{\n";
    bool first_node = true;
    xag.foreach_node([&](auto n) {
        if (xag.is_constant(n) || xag.is_pi(n)) return;

        bool has_real_candidate = false;
        for (uint32_t cidx : candidates[n]) {
            CutMatch const& c = cut_matches[cidx];
            if (c.outputs.size() > 1) continue;
            if (c.is_self_cut()) continue;
            has_real_candidate = true;
            break;
        }
        if (!has_real_candidate) return;

        if (!first_node) out << ",\n";
        first_node = false;

        bool is_and = xag.is_and(n);
        CutMatch const* sel = optimizer.selected(n);
        uint32_t sel_idx = 0;
        for (uint32_t cidx : candidates[n]) {
            CutMatch const& c = cut_matches[cidx];
            if (&c == sel && !c.is_self_cut() && c.outputs.size() == 1) {
                sel_idx = cidx;
                break;
            }
        }

        out << "    \"" << n << "\":{"
            << "\"is_and\":" << (is_and ? "true" : "false")
            << ",\"selected_cut_idx\":" << sel_idx
            << ",\"candidates\":[";
        bool first_c = true;
        for (uint32_t cidx : candidates[n]) {
            CutMatch const& c = cut_matches[cidx];
            if (c.outputs.size() > 1) continue;
            if (c.is_self_cut()) continue;
            if (!first_c) out << ',';
            first_c = false;
            uint32_t and_terms = 0;
            uint32_t xor_terms = 0;
            uint32_t max_degree = 0;
            for (auto const& a : c.anfs) {
                and_terms += static_cast<uint32_t>(a.degree2.size());
                xor_terms += static_cast<uint32_t>(a.degree1.size());
                max_degree = std::max(max_degree, a.get_degree());
            }
            out << "{\"cut_idx\":" << cidx
                << ",\"children\":"; emit_int_array(out, c.children);
            out << ",\"cost_t\":" << c.c.t_count
                << ",\"and_terms\":" << and_terms
                << ",\"xor_terms\":" << xor_terms
                << ",\"max_degree\":" << max_degree
                << ",\"and_pairs\":[";
            bool first_pair = true;
            for (auto const& a : c.anfs) {
                for (auto const& p : a.degree2) {
                    if (!first_pair) out << ',';
                    first_pair = false;
                    out << "[" << p.first << "," << p.second << "]";
                }
            }
            out << "]}";
        }
        out << "]}";
    });
    out << "\n  },\n";

    out << "  \"daomap_baseline_t\":" << run_stats.mapper.baseline_t_count << ",\n";
    out << "  \"daomap_mapped_t\":" << run_stats.mapper.mapped_t_count << "\n";
    out << "}\n";
}

}  // namespace exact_t
