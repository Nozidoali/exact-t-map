#include "network/xag.hpp"

#include <algorithm>
#include <cassert>

namespace exact_t {

XagNetwork::XagNetwork() {
    nodes_.push_back({Signal(0, false), Signal(0, false)});
}

Signal XagNetwork::get_constant(bool value) const {
    return Signal(0, value);
}

Signal XagNetwork::create_pi() {
    uint32_t idx = static_cast<uint32_t>(nodes_.size());
    nodes_.push_back({Signal(0, false), Signal(0, false)});
    num_inputs_++;
    return Signal(idx, false);
}

Signal XagNetwork::create_and(Signal a, Signal b) {
    if (a.index() > b.index()) std::swap(a, b);
    uint32_t idx = static_cast<uint32_t>(nodes_.size());
    nodes_.push_back({a, b});
    return Signal(idx, false);
}

Signal XagNetwork::create_xor(Signal a, Signal b) {
    bool output_compl = a.is_complemented() ^ b.is_complemented();
    Signal na(a.index(), false);
    Signal nb(b.index(), false);
    if (na.index() <= nb.index()) std::swap(na, nb);
    assert(na.index() != nb.index());
    uint32_t idx = static_cast<uint32_t>(nodes_.size());
    nodes_.push_back({na, nb});
    return Signal(idx, output_compl);
}

void XagNetwork::create_po(Signal s) {
    outputs_.push_back(s);
}

uint32_t XagNetwork::size() const {
    return static_cast<uint32_t>(nodes_.size());
}

uint32_t XagNetwork::num_pis() const { return num_inputs_; }

uint32_t XagNetwork::num_pos() const {
    return static_cast<uint32_t>(outputs_.size());
}

uint32_t XagNetwork::num_gates() const {
    return static_cast<uint32_t>(nodes_.size()) - 1u - num_inputs_;
}

bool XagNetwork::is_constant(uint32_t n) const { return n == 0; }

bool XagNetwork::is_pi(uint32_t n) const {
    return n >= 1 && n <= num_inputs_;
}

bool XagNetwork::is_and(uint32_t n) const {
    if (is_constant(n) || is_pi(n)) return false;
    return nodes_[n][0].index() <= nodes_[n][1].index();
}

bool XagNetwork::is_xor(uint32_t n) const {
    if (is_constant(n) || is_pi(n)) return false;
    return nodes_[n][0].index() > nodes_[n][1].index();
}

bool XagNetwork::is_complemented(Signal s) { return s.is_complemented(); }
uint32_t XagNetwork::get_node(Signal s) { return s.index(); }

TopoView::TopoView(XagNetwork const& xag) : xag_(xag) {}

}  // namespace exact_t
