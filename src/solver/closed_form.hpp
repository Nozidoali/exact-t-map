#pragma once

/*! \file closed_form.hpp
 *  \brief Closed-form T-count predictors for k<=3 bilinear cuts.
 *
 *  For k=1 single-output cuts, T_opt = 4*MC(Q) + 3 where MC(Q) is the
 *  multiplicative complexity (rank/2 of the GF(2) alternating matrix
 *  representing Q's bilinear pairs). Verified on 1792 EPFL n=7 instances
 *  by closed_form_n7.py.
 *
 *  Returns UINT32_MAX if the closed form does not apply (out of scope).
 */

#include "core/bilinear.hpp"

#include <cstdint>
#include <utility>
#include <vector>

namespace exact_t {

/*! \brief Rank of A + A^T over GF(2) where A encodes bilinear pairs. */
uint32_t gf2_alternating_rank(std::vector<std::pair<uint32_t, uint32_t>> const& pairs,
                               uint32_t n_vars);

/*! \brief MC(Q) = rank/2 for a quadratic Q on n_vars variables. */
inline uint32_t mc_quadratic(std::vector<std::pair<uint32_t, uint32_t>> const& pairs,
                              uint32_t n_vars) {
    return gf2_alternating_rank(pairs, n_vars) / 2;
}

/*! \brief T-count of a single-output bilinear cut f = y * Q(x).
 *  Returns UINT32_MAX if `bf` has more than one output or has degree-1 (k!=0) terms.
 */
uint32_t closed_form_t_single_output(BilinearFunction const& bf);

}  // namespace exact_t
