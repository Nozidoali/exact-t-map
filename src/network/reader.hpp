/*! \file reader.hpp
 *  \brief Read XAG networks from Verilog and AIGER files.
 */

#pragma once

#include "network/xag.hpp"

#include <string>

namespace exact_t {

/*! \brief Read an XAG from file, auto-detecting format by extension.
 *
 *  Supported: .v / .verilog, .aig (binary AIGER), .aag (ASCII AIGER).
 */
XagNetwork read_network(std::string const& path);

}  // namespace exact_t
