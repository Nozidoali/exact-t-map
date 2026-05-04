/*! \file mapper-stats.hpp
 *  \brief Mapping-run statistics type shared between Mapper and RunStats.
 */

#pragma once

#include <cstdint>

namespace exact_t {

/*! \brief Mapping statistics. */
struct MapperStats {
    uint32_t baseline_and_gates{0};
    uint32_t baseline_t_count{0};
    uint32_t mapped_t_count{0};
    uint32_t mapped_cnot_count{0};
    uint32_t mapped_qubit_count{0};
    double runtime_ms{0.0};
};

}  // namespace exact_t
