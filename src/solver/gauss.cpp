#include "solver/gauss.hpp"

#include <cassert>

namespace exact_t {

std::vector<int> GaussianSolver::solve(Tensor const& target, Transform const& transform) {
    std::vector<int> solution;
    uint32_t nbits = target.data().size();

    for (uint32_t r = 0; r < transform.rank; ++r) {
        Bits const& row_mask = transform.matrix[r];
        uint32_t popcount = 0;
        for (uint32_t b = 0; b < nbits; ++b) {
            if (row_mask.get(b) && target.data().get(b))
                ++popcount;
        }
        if (popcount % 2 == 1)
            solution.push_back(transform.row_to_cand[r]);
    }
    return solution;
}

std::vector<Bits> GaussianSolver::dc_freedom(Tensor const& target,
                                              Transform const& transform, uint32_t m) {
    if (!target.has_dont_cares()) return {};

    uint32_t nbits = target.data().size();
    std::vector<Bits> freedom_vecs;

    for (uint32_t dc_bit = 0; dc_bit < nbits; ++dc_bit) {
        if (target.is_care(dc_bit)) continue;
        Bits flip_vec(m);
        for (uint32_t r = 0; r < transform.rank; ++r) {
            if (transform.matrix[r].get(dc_bit))
                flip_vec.set(static_cast<uint32_t>(transform.row_to_cand[r]));
        }
        if (flip_vec.count() > 0)
            freedom_vecs.push_back(flip_vec);
    }
    return freedom_vecs;
}

void GaussianSolver::refine_basis(std::vector<Bits> const& basis, std::vector<int>& sol,
                                   uint32_t m, uint32_t max_depth) {
    if (basis.empty() || sol.empty()) return;

    Bits sol_vec(m);
    for (int idx : sol)
        sol_vec.set(static_cast<uint32_t>(idx));

    uint32_t b = static_cast<uint32_t>(basis.size());
    bool improved = true;

    while (improved) {
        improved = false;
        uint32_t best_pc = sol_vec.count();
        Bits best_flip(m);

        auto try_flip = [&](Bits const& flip) {
            Bits test = sol_vec ^ flip;
            uint32_t pc = test.count();
            if (pc < best_pc) {
                best_pc = pc;
                best_flip = flip;
                improved = true;
            }
        };

        for (uint32_t i = 0; i < b; ++i)
            try_flip(basis[i]);

        if (max_depth >= 2) {
            for (uint32_t i = 0; i < b; ++i)
                for (uint32_t j = i + 1; j < b; ++j)
                    try_flip(basis[i] ^ basis[j]);
        }

        if (max_depth >= 3) {
            for (uint32_t i = 0; i < b; ++i)
                for (uint32_t j = i + 1; j < b; ++j) {
                    Bits pair_ij = basis[i] ^ basis[j];
                    for (uint32_t k = j + 1; k < b; ++k)
                        try_flip(pair_ij ^ basis[k]);
                }
        }

        if (max_depth >= 4) {
            for (uint32_t i = 0; i < b; ++i)
                for (uint32_t j = i + 1; j < b; ++j) {
                    Bits pair_ij = basis[i] ^ basis[j];
                    for (uint32_t k = j + 1; k < b; ++k) {
                        Bits triple = pair_ij ^ basis[k];
                        for (uint32_t l = k + 1; l < b; ++l)
                            try_flip(triple ^ basis[l]);
                    }
                }
        }

        if (improved)
            sol_vec ^= best_flip;
    }

    sol.clear();
    for (uint32_t c = 0; c < m; ++c) {
        if (sol_vec.get(c))
            sol.push_back(static_cast<int>(c));
    }
}

bool GaussianSolver::verify(Phase const& phase, Tensor const& target) {
    Tensor result = Tensor::from_phase(phase);
    if (!target.has_dont_cares()) return result == target;
    uint32_t nbits = target.data().size();
    for (uint32_t b = 0; b < nbits; ++b) {
        if (!target.is_care(b)) continue;
        if (result.data().get(b) != target.data().get(b)) return false;
    }
    return true;
}

Phase GaussianSolver::to_phase(std::vector<int> const& sol, uint32_t n) {
    Phase p(n);
    for (int idx : sol)
        p.toggle_term(static_cast<uint32_t>(idx));
    return p;
}

}  // namespace exact_t
