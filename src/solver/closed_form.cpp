#include "solver/closed_form.hpp"

#include <algorithm>
#include <climits>

namespace exact_t {

uint32_t gf2_alternating_rank(std::vector<std::pair<uint32_t, uint32_t>> const& pairs,
                               uint32_t n_vars) {
    std::vector<uint64_t> rows(n_vars, 0);
    for (auto const& p : pairs) {
        uint32_t a = std::min(p.first, p.second);
        uint32_t b = std::max(p.first, p.second);
        rows[a] ^= (1ULL << b);
        rows[b] ^= (1ULL << a);
    }
    uint32_t rank = 0;
    for (uint32_t col = 0; col < n_vars; ++col) {
        uint32_t piv = n_vars;
        for (uint32_t row = rank; row < n_vars; ++row) {
            if ((rows[row] >> col) & 1ULL) { piv = row; break; }
        }
        if (piv == n_vars) continue;
        if (piv != rank) std::swap(rows[rank], rows[piv]);
        for (uint32_t row = 0; row < n_vars; ++row) {
            if (row != rank && ((rows[row] >> col) & 1ULL))
                rows[row] ^= rows[rank];
        }
        ++rank;
    }
    return rank;
}

uint32_t closed_form_t_single_output(BilinearFunction const& bf) {
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    pairs.reserve(bf.triplets().size());
    for (auto const& t : bf.triplets()) {
        if (t.k != 0) return UINT32_MAX;
        pairs.emplace_back(t.i, t.j);
    }
    uint32_t mc = mc_quadratic(pairs, bf.num_inputs());
    if (mc == 0) return 0;
    return 4 * mc + 3;
}

}  // namespace exact_t
