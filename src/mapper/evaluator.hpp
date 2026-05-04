/*! \file evaluator.hpp
 *  \brief Cost evaluator for cut matches using Clifford+T synthesis.
 */

#pragma once

#include <cstdint>

#include "mapper/bilinear.hpp"
#include "mapper/cut.hpp"
#include "solver/solver.hpp"
#include "core/bilinear.hpp"
#include "core/tensor.hpp"

#include <functional>
#include <unordered_map>

namespace exact_t {

/*! \brief Evaluates T-count cost of cut matches via synthesis. */
class CostEvaluator {
public:
    CostEvaluator(CliffordTSolver const& solver, uint32_t max_solver_n,
                  bool use_clean_ancilla, bool toffoli_mapping);

    /*! \brief Compute cost for a cut match. */
    Cost evaluate(CutMatch const& match) const;

private:
    CliffordTSolver const& solver_;
    uint32_t max_solver_n_;
    bool use_clean_ancilla_;
    bool toffoli_mapping_;

    struct CutCostKey {
        Tensor target;
        uint32_t num_outputs;
        bool operator==(CutCostKey const& o) const {
            return num_outputs == o.num_outputs && target == o.target;
        }
    };
    struct CutCostKeyHash {
        size_t operator()(CutCostKey const& k) const {
            size_t h = std::hash<Tensor>{}(k.target);
            h ^= std::hash<uint32_t>{}(k.num_outputs) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };
    mutable std::unordered_map<CutCostKey, Cost, CutCostKeyHash> cache_;
};

}  // namespace exact_t
