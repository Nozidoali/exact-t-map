#include "mapper/evaluator.hpp"

#include "solver/closed_form.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <set>

namespace exact_t {

CostEvaluator::CostEvaluator(CliffordTSolver const& solver, uint32_t max_solver_n,
                              bool use_clean_ancilla, bool toffoli_mapping)
    : solver_(solver), max_solver_n_(max_solver_n),
      use_clean_ancilla_(use_clean_ancilla), toffoli_mapping_(toffoli_mapping) {}

Cost CostEvaluator::evaluate(CutMatch const& match) const {
    if (match.degree() <= 1) return {0, 0, 0, 0, 0};

    std::set<xag_node> uo(match.outputs.begin(), match.outputs.end());
    uint32_t ni = static_cast<uint32_t>(match.children.size());
    uint32_t no = static_cast<uint32_t>(uo.size());

    if (ni + no > max_solver_n_)
        return {UINT32_MAX, UINT32_MAX, no, UINT32_MAX, UINT32_MAX};

    if (toffoli_mapping_) {
        uint32_t tc = 0;
        size_t limit = std::min(match.anfs.size(), match.outputs.size());
        for (size_t i = 0; i < limit; ++i)
            tc += static_cast<uint32_t>(match.anfs[i].degree2.size());
        return {tc, 0, no, 0, 0};
    }

    BilinearFunction bf = build_bilinear(match, match.children, uo).bf;
    uint32_t cut_n = ni + no;
    Tensor key_tensor = use_clean_ancilla_
        ? to_tensor_with_ancilla_dc(bf, solver_.n())
        : bf.to_tensor(cut_n);
    CutCostKey key{key_tensor, no};
    auto it = cache_.find(key);
    if (it != cache_.end()) return it->second;

    static bool const no_fast = std::getenv("EXACT_T_NO_CLOSED_FORM") != nullptr;
    if (!no_fast && no == 1 && !use_clean_ancilla_) {
        uint32_t closed = closed_form_t_single_output(bf);
        if (closed != UINT32_MAX) {
            uint32_t mc = (closed >= 3) ? (closed - 3) / 4 : 0;
            Cost c = {closed, 4 * mc, no, 0, 0};
            cache_.emplace(std::move(key), c);
            return c;
        }
    }

    Circuit circuit = solver_.synthesize(bf, use_clean_ancilla_);
    auto gc = circuit.count_gates();

    static char const* const diag_path = std::getenv("EXACT_T_CLOSED_FORM_LOG");
    if (diag_path != nullptr && no == 1 && !use_clean_ancilla_) {
        uint32_t closed = closed_form_t_single_output(bf);
        static FILE* diag_fp = std::fopen(diag_path, "a");
        if (diag_fp != nullptr && closed != UINT32_MAX) {
            std::fprintf(diag_fp, "n_in=%u solver_t=%u solver_cx=%u closed_t=%u gap=%d\n",
                         ni, gc.t, gc.cx, closed, (int)gc.t - (int)closed);
            std::fflush(diag_fp);
        }
    }

    Cost c = {gc.t, gc.cx, no, gc.s, gc.z};
    cache_.emplace(std::move(key), c);
    return c;
}

}  // namespace exact_t
