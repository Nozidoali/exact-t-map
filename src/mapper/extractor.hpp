/*! \file extractor.hpp
 *  \brief Circuit extraction from optimized DAOMap selections.
 */

#pragma once

#include <cstdint>

#include "mapper/allocate.hpp"
#include "mapper/daomap.hpp"
#include "core/bilinear.hpp"
#include "solver/solver.hpp"

#include <map>
#include <set>
#include <vector>

namespace exact_t {

/*! \brief Extracts Clifford+T circuits from DAOMap-optimized XAG. */
class CircuitExtractor {
  public:
    CircuitExtractor(xag_network const& xag, DaoMapOptimizer const& optimizer,
                      CliffordTSolver const& solver, bool toffoli_mapping,
                      bool use_clean_ancilla = false);

    /*! \brief Extract full circuit with given allocation strategy. */
    Circuit extract(AllocationStrategy allocation);

  private:
    Circuit extract_circuit_recursive(xag_node node, QubitAllocator& alloc);
    size_t find_anf_index(CutMatch const& match, xag_node target) const;
    std::vector<xag_node> extract_children(CutMatch const& match,
                                            QubitAllocator& alloc, Circuit& result);
    void extract_toffoli_circuit(CutMatch const& match, xag_node node,
                                  ANF2<xag_node> const& anf,
                                  std::set<xag_node> const& unique_outputs,
                                  QubitAllocator& alloc, Circuit& result);
    void extract_solver_circuit(CutMatch const& match, xag_node node,
                                 std::vector<xag_node> const& input_nodes,
                                 std::set<xag_node> const& unique_outputs,
                                 QubitAllocator& alloc, Circuit& result);
    Circuit synthesize_degree_zero(xag_node node, ANF2<xag_node> const& anf,
                                    QubitAllocator& alloc);
    Circuit synthesize_degree_one(xag_node node, CutMatch const& match,
                                   ANF2<xag_node> const& anf, QubitAllocator& alloc);

    xag_network const& xag_;
    DaoMapOptimizer const& optimizer_;
    CliffordTSolver const& solver_;
    bool toffoli_mapping_;
    bool use_clean_ancilla_;
};

}  // namespace exact_t
