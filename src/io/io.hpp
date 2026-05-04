/*! \file io.hpp
 *  \brief Read XAG networks from various file formats.
 */

#pragma once

#include "network/xag.hpp"

#include <string>

namespace exact_t {

/*! \brief Read an XAG network from file with auto-detection by extension.
 *
 *  Supported formats: .v/.verilog, .aig, .aag
 */
XagNetwork read_network(std::string const& path);

}  // namespace exact_t
