/*! \file allocate.hpp
 *  \brief Qubit allocator for XAG to quantum circuit extraction.
 */

#pragma once

#include "mapper/cut.hpp"

#include <cassert>
#include <cstdint>
#include <functional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace exact_t {

/*! \brief Qubit allocation strategy for XOR nodes. */
enum class AllocationStrategy { NAIVE, BASIC, LTFI };

/*! \brief Tracks a borrowed qubit for deferred uncompute. */
struct BorrowInfo {
    xag_node xor_node;
    xag_node borrowed_child;
    uint32_t qubit;
    ANF2<xag_node> anf;
};

/*! \brief Qubit allocator for XAG extraction. */
class QubitAllocator {
  public:
    QubitAllocator(xag_network const& xag,
                   std::function<CutMatch const*(xag_node)> get_selected_cut,
                   AllocationStrategy strategy = AllocationStrategy::LTFI);

    virtual ~QubitAllocator() = default;

    uint32_t assign_pi(xag_node node);
    virtual void assign_node(xag_node node, CutMatch const& match);
    void assign_outputs(xag_node node, CutMatch const& match,
                        std::set<xag_node>& unique_outputs);
    uint32_t allocate_ancilla();
    uint32_t ensure_qubit(xag_node node);
    uint32_t get_qubit(xag_node node) const;
    bool has_qubit(xag_node node) const;
    void mark_visited(xag_node node);
    bool is_visited(xag_node node) const;
    void consume_refs(CutMatch const& match);
    bool is_primary_output(xag_node node) const;
    uint32_t next_qubit() const;
    uint32_t qubits_reused() const;
    uint32_t qubits_forced() const;
    uint32_t qubits_borrowed() const;
    void record_borrow_anf(xag_node node, ANF2<xag_node> const& anf);
    std::vector<BorrowInfo> get_pending_uncomputes();

  protected:
    virtual bool can_reuse(xag_node node, CutMatch const& match,
                           uint32_t& reuse_qubit);
    bool can_borrow(xag_node node, CutMatch const& match,
                    uint32_t& borrow_qubit, xag_node& borrowed_child);

    xag_network const& xag_;
    std::function<CutMatch const*(xag_node)> get_selected_cut_;
    std::unordered_map<xag_node, uint32_t> node_to_qubit_;
    uint32_t next_qubit_idx_{0};
    std::unordered_set<xag_node> visited_;
    uint32_t qubits_reused_{0};
    uint32_t qubits_forced_{0};
    uint32_t qubits_borrowed_{0};
    AllocationStrategy strategy_;
    std::unordered_set<xag_node> primary_outputs_;
    std::unordered_map<xag_node, uint32_t> remaining_refs_;
    std::unordered_map<xag_node, BorrowInfo> active_borrows_;
};

}  // namespace exact_t
