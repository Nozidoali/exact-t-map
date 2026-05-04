#include "mapper/extractor.hpp"
#include "mapper/bilinear.hpp"

#include <algorithm>
#include <set>

namespace exact_t {

namespace {

void append_identity(Circuit& dest, Circuit const& src) {
    if (src.gates().empty()) return;
    uint32_t nq = src.num_qubits();
    for (auto const& g : src.gates()) {
        nq = std::max(nq, g.target + 1);
        if (g.type == GateType::CX) nq = std::max(nq, g.control + 1);
    }
    std::vector<uint32_t> mapping(nq);
    for (uint32_t i = 0; i < nq; ++i) mapping[i] = i;
    dest.append(src, mapping);
    dest.set_num_qubits(std::max(dest.num_qubits(), nq));
}

}  // namespace

CircuitExtractor::CircuitExtractor(xag_network const& xag,
                                     DaoMapOptimizer const& optimizer,
                                     CliffordTSolver const& solver,
                                     bool toffoli_mapping,
                                     bool use_clean_ancilla)
    : xag_(xag), optimizer_(optimizer), solver_(solver),
      toffoli_mapping_(toffoli_mapping), use_clean_ancilla_(use_clean_ancilla) {}

Circuit CircuitExtractor::extract(AllocationStrategy allocation) {
    Circuit result;
    QubitAllocator alloc(xag_, [this](xag_node n) { return optimizer_.selected(n); },
                         allocation);

    xag_.foreach_pi([&](auto pi_node) {
        result.add_input(alloc.assign_pi(pi_node));
    });

    xag_.foreach_po([&](auto po_signal) {
        xag_node po_node = xag_.get_node(po_signal);
        bool po_complement = xag_.is_complemented(po_signal);
        Circuit po_circuit = extract_circuit_recursive(po_node, alloc);
        append_identity(result, po_circuit);
        if (alloc.has_qubit(po_node)) {
            uint32_t oq = alloc.get_qubit(po_node);
            result.add_output(oq);
            if (po_complement) result.x(oq);
        }
    });

    result.finalize_qubit_count();
    return result;
}

size_t CircuitExtractor::find_anf_index(CutMatch const& match,
                                          xag_node target) const {
    for (size_t i = 0; i < match.outputs.size(); ++i)
        if (match.outputs[i] == target) return i;
    return 0;
}

std::vector<xag_node> CircuitExtractor::extract_children(
    CutMatch const& match, QubitAllocator& alloc, Circuit& result) {
    std::vector<xag_node> input_nodes;
    for (auto child_node : match.children) {
        Circuit child = extract_circuit_recursive(child_node, alloc);
        append_identity(result, child);
        input_nodes.push_back(child_node);
        alloc.ensure_qubit(child_node);
    }
    return input_nodes;
}

void CircuitExtractor::extract_toffoli_circuit(
    CutMatch const& match, xag_node node, ANF2<xag_node> const& anf,
    std::set<xag_node> const& unique_outputs,
    QubitAllocator& alloc, Circuit& result) {
    for (auto out : unique_outputs) alloc.mark_visited(out);

    auto emit = [&](xag_node out_node, ANF2<xag_node> const& a) {
        if (!alloc.has_qubit(out_node)) return;
        uint32_t out_q = alloc.get_qubit(out_node);
        for (auto const& p : a.degree2) {
            if (!alloc.has_qubit(p.first) || !alloc.has_qubit(p.second)) continue;
            result.toffoli(alloc.get_qubit(p.first), alloc.get_qubit(p.second), out_q);
        }
        for (auto t : a.degree1) {
            if (!alloc.has_qubit(t)) continue;
            result.cx(alloc.get_qubit(t), out_q);
        }
        if (a.degree0) result.x(out_q);
    };

    if (!match.outputs.empty() && !match.anfs.empty()) {
        size_t limit = std::min(match.outputs.size(), match.anfs.size());
        for (size_t i = 0; i < limit; ++i)
            emit(match.outputs[i], match.anfs[i]);
    } else {
        emit(node, anf);
    }
    result.set_num_qubits(std::max(result.num_qubits(), alloc.next_qubit()));
}

void CircuitExtractor::extract_solver_circuit(
    CutMatch const& match, xag_node node,
    std::vector<xag_node> const& input_nodes,
    std::set<xag_node> const& unique_outputs,
    QubitAllocator& alloc, Circuit& result) {
    auto [bf, node_to_bf_input, node_to_bf_output] =
        build_bilinear(match, input_nodes, unique_outputs);

    std::map<xag_node, std::set<uint32_t>> degree1_terms;
    std::set<xag_node> outputs_to_negate;

    for (size_t i = 0; i < match.anfs.size() && i < match.outputs.size(); ++i) {
        xag_node t_output = match.outputs[i];
        for (auto term : match.anfs[i].degree1)
            if (node_to_bf_input.count(term) && node_to_bf_output.count(t_output))
                degree1_terms[term].insert(node_to_bf_output[t_output]);
        if (match.anfs[i].degree0) outputs_to_negate.insert(t_output);
    }

    Circuit cut_circuit = solver_.synthesize(bf, use_clean_ancilla_);

    std::set<uint32_t> used_qubits;
    for (auto const& gate : cut_circuit.gates()) {
        used_qubits.insert(gate.target);
        if (gate.type == GateType::CX) used_qubits.insert(gate.control);
    }

    std::vector<uint32_t> mapping(cut_circuit.num_qubits());
    std::unordered_map<uint32_t, uint32_t> compact_map;

    for (uint32_t i = 0; i < bf.num_inputs(); ++i)
        compact_map[i] = alloc.get_qubit(input_nodes[i]);
    uint32_t oi = 0;
    for (auto out : unique_outputs)
        compact_map[bf.num_inputs() + oi++] = alloc.get_qubit(out);

    for (uint32_t i = 0; i < cut_circuit.num_qubits(); ++i) {
        if (compact_map.count(i))
            mapping[i] = compact_map[i];
        else if (used_qubits.count(i))
            mapping[i] = alloc.allocate_ancilla();
        else
            mapping[i] = 0;
    }

    result.append(cut_circuit, mapping);
    result.set_num_qubits(std::max(result.num_qubits(), alloc.next_qubit()));

    for (auto out : outputs_to_negate)
        if (alloc.has_qubit(out)) result.x(alloc.get_qubit(out));

    for (auto const& [input_node, output_bf_indices] : degree1_terms) {
        if (!alloc.has_qubit(input_node)) continue;
        uint32_t iq = alloc.get_qubit(input_node);
        for (uint32_t bf_out_idx : output_bf_indices) {
            auto it = unique_outputs.begin();
            std::advance(it, bf_out_idx);
            if (alloc.has_qubit(*it)) result.cx(iq, alloc.get_qubit(*it));
        }
    }
}

Circuit CircuitExtractor::synthesize_degree_zero(xag_node node,
                                                   ANF2<xag_node> const& anf,
                                                   QubitAllocator& alloc) {
    Circuit result;
    if (anf.degree0) result.x(alloc.get_qubit(node));
    return result;
}

Circuit CircuitExtractor::synthesize_degree_one(xag_node node,
                                                  CutMatch const& match,
                                                  ANF2<xag_node> const& anf,
                                                  QubitAllocator& alloc) {
    Circuit result;
    for (auto child_node : match.children) {
        Circuit child = extract_circuit_recursive(child_node, alloc);
        append_identity(result, child);
    }
    alloc.assign_node(node, match);
    alloc.record_borrow_anf(node, anf);

    for (auto term : anf.degree1) {
        if (alloc.has_qubit(term)) {
            uint32_t src = alloc.get_qubit(term);
            uint32_t tgt = alloc.get_qubit(node);
            if (src != tgt) result.cx(src, tgt);
        }
    }
    if (anf.degree0) result.x(alloc.get_qubit(node));
    return result;
}

Circuit CircuitExtractor::extract_circuit_recursive(xag_node node,
                                                      QubitAllocator& alloc) {
    Circuit result;
    if (alloc.is_visited(node)) return result;
    if (xag_.is_constant(node) || xag_.is_pi(node)) {
        alloc.mark_visited(node);
        alloc.ensure_qubit(node);
        return result;
    }
    alloc.mark_visited(node);

    CutMatch const* match_ptr = optimizer_.selected(node);
    if (!match_ptr) return result;
    CutMatch const& match = *match_ptr;

    size_t anf_idx = find_anf_index(match, node);
    if (anf_idx >= match.anfs.size()) anf_idx = 0;
    ANF2<xag_node> const& anf = match.anfs.empty() ? ANF2<xag_node>() : match.anfs[anf_idx];

    if (match.degree() != 1) alloc.assign_node(node, match);

    if (match.degree() == 0) {
        result = synthesize_degree_zero(node, anf, alloc);
    } else if (match.degree() == 1) {
        result = synthesize_degree_one(node, match, anf, alloc);
    } else if (anf.degree2.empty()) {
        for (auto term : anf.degree1)
            if (alloc.has_qubit(term))
                result.cx(alloc.get_qubit(term), alloc.get_qubit(node));
        if (anf.degree0) result.x(alloc.get_qubit(node));
    } else {
        auto input_nodes = extract_children(match, alloc, result);
        std::set<xag_node> unique_outputs;
        alloc.assign_outputs(node, match, unique_outputs);
        if (!unique_outputs.empty()) {
            if (toffoli_mapping_)
                extract_toffoli_circuit(match, node, anf, unique_outputs, alloc, result);
            else
                extract_solver_circuit(match, node, input_nodes, unique_outputs, alloc, result);
        }
        if (match.outputs.size() > 1)
            for (auto out : match.outputs) alloc.mark_visited(out);
    }

    alloc.consume_refs(match);
    for (auto const& info : alloc.get_pending_uncomputes()) {
        if (info.anf.degree0) result.x(info.qubit);
        for (auto const& term : info.anf.degree1)
            if (alloc.has_qubit(term)) {
                uint32_t src = alloc.get_qubit(term);
                if (src != info.qubit) result.cx(src, info.qubit);
            }
    }
    return result;
}

}  // namespace exact_t
