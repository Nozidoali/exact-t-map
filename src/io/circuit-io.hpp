/*! \file circuit-io.hpp
 *  \brief QASM and QC serialization for Clifford+T circuits.
 */

#pragma once

#include "core/circuit.hpp"
#include <string>

namespace exact_t {

/*! \brief Serialize circuit to OpenQASM 2.0 string. */
std::string to_qasm(Circuit const& circuit);

/*! \brief Serialize circuit to QC format string. */
std::string to_qc(Circuit const& circuit);

/*! \brief Write circuit to file in QASM format. */
void write_qasm(Circuit const& circuit, std::string const& path);

/*! \brief Write circuit to file in QC format. */
void write_qc(Circuit const& circuit, std::string const& path);

}  // namespace exact_t
