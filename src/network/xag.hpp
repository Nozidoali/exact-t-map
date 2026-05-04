/*! \file xag.hpp
 *  \brief Self-contained XAG network, Signal, NodeMap, and TopoView.
 *
 *  Provides the same API surface as mockturtle's xag_network, node_map,
 *  and topo_view, restricted to the methods used by exact-t.
 */

#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace exact_t {

/*! \brief Signal in an XAG network.
 *
 *  Bit 0 encodes complement, bits 1..31 encode the node index.
 */
class Signal {
    uint32_t data_{0};

  public:
    Signal() = default;
    Signal(uint32_t index, bool complement)
        : data_((index << 1) | static_cast<uint32_t>(complement)) {}

    uint32_t index() const { return data_ >> 1; }
    bool is_complemented() const { return data_ & 1u; }

    Signal operator!() const {
        Signal s;
        s.data_ = data_ ^ 1u;
        return s;
    }

    bool operator==(Signal o) const { return data_ == o.data_; }
    bool operator!=(Signal o) const { return data_ != o.data_; }
    bool operator<(Signal o) const { return data_ < o.data_; }
};

/*! \brief XOR-AND graph network.
 *
 *  Node 0 is the constant-false node.
 *  Nodes 1 .. num_inputs_ are primary inputs.
 *  Remaining nodes are gates.
 *
 *  Gate type convention (matching mockturtle):
 *    AND  — fanin[0].index() <= fanin[1].index()
 *    XOR  — fanin[0].index() >  fanin[1].index()
 */
class XagNetwork {
  public:
    using signal = Signal;
    using node = uint32_t;

    XagNetwork();

    Signal get_constant(bool value) const;
    Signal create_pi();
    Signal create_and(Signal a, Signal b);
    Signal create_xor(Signal a, Signal b);
    void create_po(Signal s);

    uint32_t size() const;
    uint32_t num_pis() const;
    uint32_t num_pos() const;
    uint32_t num_gates() const;

    bool is_constant(uint32_t n) const;
    bool is_pi(uint32_t n) const;
    bool is_and(uint32_t n) const;
    bool is_xor(uint32_t n) const;

    static bool is_complemented(Signal s);
    static uint32_t get_node(Signal s);

    template <typename Fn>
    void foreach_node(Fn&& fn) const {
        for (uint32_t i = 0; i < static_cast<uint32_t>(nodes_.size()); ++i)
            fn(i);
    }

    template <typename Fn>
    void foreach_pi(Fn&& fn) const {
        for (uint32_t i = 1; i <= num_inputs_; ++i)
            fn(i);
    }

    template <typename Fn>
    void foreach_po(Fn&& fn) const {
        for (auto const& s : outputs_)
            fn(s);
    }

    template <typename Fn>
    void foreach_fanin(uint32_t n, Fn&& fn) const {
        if (!is_constant(n) && !is_pi(n)) {
            fn(nodes_[n][0]);
            fn(nodes_[n][1]);
        }
    }

  private:
    std::vector<std::array<Signal, 2>> nodes_;
    uint32_t num_inputs_{0};
    std::vector<Signal> outputs_;
};

/*! \brief Map from node IDs to values, indexed by uint32_t. */
template <typename T>
class NodeMap {
    std::vector<T> data_;

  public:
    explicit NodeMap(XagNetwork const& xag) : data_(xag.size()) {}

    T& operator[](uint32_t n) { return data_[n]; }
    T const& operator[](uint32_t n) const { return data_[n]; }
};

/*! \brief Topological-order view over an XAG network.
 *
 *  Since nodes are created in dependency order (fanins always have lower
 *  indices), the natural index order 0 .. size()-1 is already topological.
 */
class TopoView {
    XagNetwork const& xag_;

  public:
    explicit TopoView(XagNetwork const& xag);

    template <typename Fn>
    void foreach_node(Fn&& fn) const {
        for (uint32_t i = 0; i < xag_.size(); ++i)
            fn(i);
    }
};

using xag_network = XagNetwork;
using xag_signal = Signal;
using xag_node = uint32_t;
template <typename T> using node_map = NodeMap<T>;
using topo_view = TopoView;

}  // namespace exact_t
