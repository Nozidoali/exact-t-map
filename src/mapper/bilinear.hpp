/*! \file bilinear.hpp
 *  \brief Shared BilinearBuild struct and builder for mapper.
 */

#pragma once

#include <cstdint>

#include "mapper/cut.hpp"
#include "core/bilinear.hpp"

#include <set>
#include <unordered_map>
#include <vector>

namespace exact_t {

/*! \brief Result of building a BilinearFunction from a cut match. */
struct BilinearBuild {
    BilinearFunction bf;
    std::unordered_map<xag_node, uint32_t> node_to_input;
    std::unordered_map<xag_node, uint32_t> node_to_output;
};

/*! \brief Construct a BilinearFunction from a cut match. */
BilinearBuild build_bilinear(CutMatch const& match,
                              std::vector<xag_node> const& input_nodes,
                              std::set<xag_node> const& unique_outputs);

}  // namespace exact_t
