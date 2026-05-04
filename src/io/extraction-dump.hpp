/*! \file extraction-dump.hpp
 *  \brief Serialize cut candidate set + current DaoMap selection to JSON
 *         for offline ILP-based extraction analysis.
 */

#pragma once

#include <string>

namespace exact_t {

class Mapper;
class DaoMapOptimizer;

/*! \brief Dump the extraction problem (POs, non-PI nodes, their single-output
 *         candidate cuts, cut children, cut costs, and DaoMap's current
 *         selection) to a JSON file.
 *
 *  Cuts with `outputs.size() > 1` are filtered out (Phase 1 skips multi-output).
 *  XOR-node cuts are included with `cost_t = 0` so the ILP can correctly enforce
 *  coverage constraints on their children.
 *
 *  Asserts on file open failure (consistent with other io/ writers).
 */
void write_extraction_dump(Mapper const& mapper, DaoMapOptimizer const& optimizer,
                           std::string const& path);

}  // namespace exact_t
