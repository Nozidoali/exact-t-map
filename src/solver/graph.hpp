/*! \file graph.hpp
 *  \brief Graph cache for parity synthesis paths.
 */

#pragma once

#include "core/circuit.hpp"
#include "solver/parity-matrix.hpp"
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>

namespace exact_t {

/*! \brief Edge in parity synthesis graph. */
struct Edge {
    uint32_t next_id{UINT32_MAX};
    bool is_phase_gate{false};
    uint32_t control{0};
    uint32_t target{0};
    uint32_t phase_row{0};
    uint32_t phase_col{0};
    std::vector<uint32_t> row_perm;
    std::vector<uint32_t> col_perm;

    /*! \brief Apply edge transformation to matrix. */
    ParityMatrix operator()(ParityMatrix const& pmat) const;
};

/*! \brief Node in parity synthesis graph. */
struct GraphNode {
    uint32_t id{0};
    std::shared_ptr<ParityMatrix> matrix;
    Edge edge;
    uint32_t cnot_cost{UINT32_MAX};

    GraphNode() = default;
    GraphNode(uint32_t nid, ParityMatrix const& m);
};

/*! \brief Graph cache for parity synthesis paths. */
class Graph {
  public:
    Graph() = default;
    explicit Graph(uint32_t n) : n_(n) {}

    uint32_t get_n() const { return n_; }
    void set_n(uint32_t n) { n_ = n; }

    /*! \brief Add node for matrix state. */
    uint32_t add_node(ParityMatrix const& m);

    /*! \brief Find existing node for matrix. */
    uint32_t find_node(ParityMatrix const& m) const;

    /*! \brief Get node by ID. */
    GraphNode const* get_node(uint32_t id) const;

    /*! \brief Set edge between nodes. */
    void set_edge(uint32_t src, uint32_t tgt, bool is_phase,
                  uint32_t ctrl, uint32_t tgt_q,
                  uint32_t phase_row, uint32_t phase_col,
                  std::vector<uint32_t> const& row_perm,
                  std::vector<uint32_t> const& col_perm);

    /*! \brief Set CNOT cost for node. */
    void set_node_cost(uint32_t id, uint32_t cnot_cost);

    size_t node_count() const { return nodes_.size(); }

    /*! \brief Extract circuit from graph path. */
    Circuit extract_circuit(uint32_t start_id, ParityMatrix const& input_pmat) const;

    /*! \brief Load graph from file. */
    void load(std::string const& path);

    /*! \brief Save graph to file. */
    void save(std::string const& path) const;

    bool is_dirty() const { return dirty_; }
    void set_clean() { dirty_ = false; }

    /*! \brief Get graph for n qubits. */
    Graph& get_graph(uint32_t n);
    Graph const& get_graph(uint32_t n) const;

    void record_lookup() const { lookups_++; }
    void record_hit() const { hits_++; }
    void record_miss() const { misses_++; }

    void load_file(std::string const& base_path);
    void save_file(std::string const& base_path);

    uint64_t get_lookups() const { return lookups_; }
    uint64_t get_hits() const { return hits_; }
    uint64_t get_misses() const { return misses_; }

    bool is_any_dirty() const;
    void set_all_clean();

  private:
    uint32_t n_{0};
    std::unordered_map<ParityMatrix, uint32_t, ParityMatrixHash> mat_to_id_;
    std::vector<GraphNode> nodes_;
    mutable bool dirty_{false};

    std::unordered_map<uint32_t, std::unique_ptr<Graph>> graphs_;
    mutable uint64_t lookups_{0};
    mutable uint64_t hits_{0};
    mutable uint64_t misses_{0};
};

}  // namespace exact_t
