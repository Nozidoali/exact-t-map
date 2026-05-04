/*! \file gauss.hpp
 *  \brief Gaussian elimination solver and basis refinement for T-count optimization.
 */

#pragma once

#include <cstdint>
#include <vector>

#include "core/bits.hpp"
#include "core/phase.hpp"
#include "core/tensor.hpp"
#include "solver/solver.hpp"

namespace exact_t {

/*! \brief Gaussian elimination solver with basis refinement. */
class GaussianSolver {
public:
    static std::vector<int> solve(Tensor const& target, Transform const& transform);
    static std::vector<Bits> dc_freedom(Tensor const& target, Transform const& transform, uint32_t m);
    static void refine_basis(std::vector<Bits> const& basis, std::vector<int>& sol, uint32_t m, uint32_t max_depth);
    static bool verify(Phase const& phase, Tensor const& target);
    static Phase to_phase(std::vector<int> const& sol, uint32_t n);
};

}  // namespace exact_t
