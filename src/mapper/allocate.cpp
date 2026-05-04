#include "mapper/allocate.hpp"

#include <queue>

namespace exact_t {

QubitAllocator::QubitAllocator(
    xag_network const& xag,
    std::function<CutMatch const*(xag_node)> get_selected_cut,
    AllocationStrategy strategy)
    : xag_(xag), get_selected_cut_(std::move(get_selected_cut)),
      strategy_(strategy) {
    xag_.foreach_po([&](auto po_signal) {
        primary_outputs_.insert(xag_.get_node(po_signal));
    });

    std::set<xag_node> reachable;
    std::queue<xag_node> worklist;
    xag_.foreach_po([&](auto po_signal) {
        xag_node po_node = xag_.get_node(po_signal);
        if (reachable.insert(po_node).second)
            worklist.push(po_node);
    });
    while (!worklist.empty()) {
        xag_node n = worklist.front();
        worklist.pop();
        if (xag_.is_constant(n) || xag_.is_pi(n)) continue;
        CutMatch const* match = get_selected_cut_(n);
        if (!match) continue;
        for (xag_node child : match->children)
            if (reachable.insert(child).second) worklist.push(child);
    }

    std::unordered_set<CutMatch const*> processed_cuts;
    for (xag_node n : reachable) {
        if (xag_.is_constant(n) || xag_.is_pi(n)) continue;
        CutMatch const* match = get_selected_cut_(n);
        if (!match) continue;
        if (match->outputs.size() > 1 && !processed_cuts.insert(match).second)
            continue;
        for (xag_node child : match->children)
            remaining_refs_[child]++;
    }
    xag_.foreach_po([&](auto po_signal) {
        remaining_refs_[xag_.get_node(po_signal)]++;
    });
}

uint32_t QubitAllocator::assign_pi(xag_node node) {
    uint32_t qubit = next_qubit_idx_++;
    node_to_qubit_[node] = qubit;
    return qubit;
}

void QubitAllocator::assign_node(xag_node node, CutMatch const& match) {
    if (node_to_qubit_.count(node)) return;
    uint32_t reuse_qubit;
    if (can_reuse(node, match, reuse_qubit)) {
        node_to_qubit_[node] = reuse_qubit;
        qubits_reused_++;
        return;
    }
    uint32_t borrow_qubit;
    xag_node borrowed_child;
    if (can_borrow(node, match, borrow_qubit, borrowed_child)) {
        node_to_qubit_[node] = borrow_qubit;
        qubits_borrowed_++;
        active_borrows_[node] = {node, borrowed_child, borrow_qubit, {}};
        return;
    }
    node_to_qubit_[node] = next_qubit_idx_++;
    if (match.degree() == 1 && strategy_ != AllocationStrategy::NAIVE)
        qubits_forced_++;
}

bool QubitAllocator::can_reuse(xag_node, CutMatch const& match,
                                uint32_t& reuse_qubit) {
    if (strategy_ == AllocationStrategy::NAIVE) return false;
    if (match.degree() != 1 || match.children.empty()) return false;

    if (strategy_ == AllocationStrategy::BASIC) {
        xag_node child = match.children[0];
        auto ref_it = remaining_refs_.find(child);
        uint32_t refs = ref_it != remaining_refs_.end() ? ref_it->second : 0;
        if (refs > 1 || !node_to_qubit_.count(child)) return false;
        if (primary_outputs_.count(child)) return false;
        reuse_qubit = node_to_qubit_.at(child);
        return true;
    }

    for (xag_node child : match.children) {
        auto ref_it = remaining_refs_.find(child);
        uint32_t refs = ref_it != remaining_refs_.end() ? ref_it->second : 0;
        if (refs <= 1 && node_to_qubit_.count(child)) {
            if (primary_outputs_.count(child)) continue;
            reuse_qubit = node_to_qubit_.at(child);
            return true;
        }
    }
    return false;
}

bool QubitAllocator::can_borrow(xag_node node, CutMatch const& match,
                                 uint32_t& borrow_qubit,
                                 xag_node& borrowed_child) {
    if (strategy_ == AllocationStrategy::NAIVE) return false;
    if (match.degree() != 1 || match.children.empty()) return false;
    if (primary_outputs_.count(node)) return false;

    for (xag_node child : match.children) {
        if (!primary_outputs_.count(child)) continue;
        if (!node_to_qubit_.count(child)) continue;
        borrow_qubit = node_to_qubit_.at(child);
        borrowed_child = child;
        return true;
    }
    return false;
}

void QubitAllocator::assign_outputs(xag_node node, CutMatch const& match,
                                     std::set<xag_node>& unique_outputs) {
    if (!match.outputs.empty())
        unique_outputs.insert(match.outputs.begin(), match.outputs.end());
    else if (primary_outputs_.count(node))
        unique_outputs.insert(node);

    for (xag_node out : unique_outputs) {
        if (!node_to_qubit_.count(out)) {
            assert(!visited_.count(out));
            node_to_qubit_[out] = next_qubit_idx_++;
        }
    }
}

uint32_t QubitAllocator::allocate_ancilla() { return next_qubit_idx_++; }

uint32_t QubitAllocator::ensure_qubit(xag_node node) {
    if (!node_to_qubit_.count(node))
        node_to_qubit_[node] = next_qubit_idx_++;
    return node_to_qubit_[node];
}

uint32_t QubitAllocator::get_qubit(xag_node node) const {
    auto it = node_to_qubit_.find(node);
    assert(it != node_to_qubit_.end());
    return it->second;
}

bool QubitAllocator::has_qubit(xag_node node) const {
    return node_to_qubit_.count(node) > 0;
}

void QubitAllocator::mark_visited(xag_node node) { visited_.insert(node); }
bool QubitAllocator::is_visited(xag_node node) const { return visited_.count(node) > 0; }

void QubitAllocator::consume_refs(CutMatch const& match) {
    for (xag_node child : match.children) {
        auto it = remaining_refs_.find(child);
        if (it != remaining_refs_.end() && it->second > 0)
            it->second--;
    }
}

bool QubitAllocator::is_primary_output(xag_node node) const {
    return primary_outputs_.count(node) > 0;
}

uint32_t QubitAllocator::next_qubit() const { return next_qubit_idx_; }
uint32_t QubitAllocator::qubits_reused() const { return qubits_reused_; }
uint32_t QubitAllocator::qubits_forced() const { return qubits_forced_; }
uint32_t QubitAllocator::qubits_borrowed() const { return qubits_borrowed_; }

void QubitAllocator::record_borrow_anf(xag_node node,
                                        ANF2<xag_node> const& anf) {
    auto it = active_borrows_.find(node);
    if (it != active_borrows_.end()) it->second.anf = anf;
}

std::vector<BorrowInfo> QubitAllocator::get_pending_uncomputes() {
    std::vector<BorrowInfo> result;
    auto it = active_borrows_.begin();
    while (it != active_borrows_.end()) {
        auto ref_it = remaining_refs_.find(it->first);
        uint32_t refs = ref_it != remaining_refs_.end() ? ref_it->second : 0;
        if (refs == 0) {
            result.push_back(it->second);
            it = active_borrows_.erase(it);
        } else {
            ++it;
        }
    }
    return result;
}

}  // namespace exact_t
