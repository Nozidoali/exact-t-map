/*! \file affine.hpp
 *  \brief Affine A*-based phase polynomial optimizer.
 */

#pragma once

#include "core/circuit.hpp"
#include "core/linear.hpp"
#include "solver/solver.hpp"

#include <cstdint>
#include <utility>
#include <vector>

namespace exact_t {
namespace affine {

/*! \brief Parameters for Affine optimization. */
struct AffineParams {
    uint32_t max_queue_size{10000};
    uint32_t max_block_qubits{12};
    uint32_t min_t_count{2};
    double timeout_s{60.0};
    bool verbose{false};
};

/*! \brief Statistics from Affine optimization. */
struct AffineStats {
    uint32_t initial_cnot_count{0};
    uint32_t final_cnot_count{0};
    uint32_t blocks_processed{0};
    uint32_t states_explored{0};
    uint32_t timeouts{0};
    double total_time_ms{0.0};
};

/*! \brief Result of parity synthesis using A* search. */
struct AffineResult {
    std::vector<std::pair<uint32_t, uint32_t>> cnots;
    std::vector<std::pair<size_t, std::pair<uint32_t, int>>> phases;
    uint32_t cnot_count{0};
    uint32_t states_explored{0};
    bool success{false};
};

/*! \brief Synthesize CNOT sequence using A* search. */
AffineResult astar_synthesize(Z8Phase const& phase, BitMatrix const& target,
                                  uint32_t n, AffineParams const& params);

/*! \brief Convert synthesis result to circuit. */
Circuit result_to_circuit(AffineResult const& result,
                           Z8Phase const& phase, uint32_t n);

/*! \brief Fallback Gray code synthesis. */
AffineResult gray_fallback(Z8Phase const& phase, BitMatrix const& target,
                               uint32_t n);

}  // namespace affine
}  // namespace exact_t
