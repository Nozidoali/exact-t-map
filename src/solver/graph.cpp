#include "solver/graph.hpp"

#include <fstream>
#include <sstream>
#include <unordered_set>

namespace exact_t {

ParityMatrix Edge::operator()(ParityMatrix const& pmat) const {
    ParityMatrix result = pmat.copy();
    if (is_phase_gate) {
        if (phase_col < result.m() && result.is_column_one_hot(phase_col) &&
            result.get_one_hot_row(phase_col) == phase_row)
            result.remove_column(phase_col);
    } else {
        if (control < result.n() && target < result.n() && control != target)
            result.apply_cnot(control, target);
    }
    std::vector<uint32_t> rp, cp;
    return result.canonicalize(rp, cp);
}

GraphNode::GraphNode(uint32_t nid, ParityMatrix const& m)
    : id(nid), matrix(std::make_shared<ParityMatrix>(m)), cnot_cost(UINT32_MAX) {}

uint32_t Graph::add_node(ParityMatrix const& m) {
    auto it = mat_to_id_.find(m);
    if (it != mat_to_id_.end()) return it->second;
    uint32_t id = static_cast<uint32_t>(nodes_.size());
    nodes_.push_back(GraphNode(id, m));
    mat_to_id_[m] = id;
    dirty_ = true;
    return id;
}

uint32_t Graph::find_node(ParityMatrix const& m) const {
    auto it = mat_to_id_.find(m);
    return it != mat_to_id_.end() ? it->second : UINT32_MAX;
}

GraphNode const* Graph::get_node(uint32_t id) const {
    return id < nodes_.size() ? &nodes_[id] : nullptr;
}

void Graph::set_edge(uint32_t src, uint32_t tgt, bool is_phase,
                     uint32_t ctrl, uint32_t tgt_q,
                     uint32_t phase_row, uint32_t phase_col,
                     std::vector<uint32_t> const& row_perm,
                     std::vector<uint32_t> const& col_perm) {
    if (src >= nodes_.size()) return;
    nodes_[src].edge.next_id = tgt;
    nodes_[src].edge.is_phase_gate = is_phase;
    nodes_[src].edge.control = ctrl;
    nodes_[src].edge.target = tgt_q;
    nodes_[src].edge.phase_row = phase_row;
    nodes_[src].edge.phase_col = phase_col;
    nodes_[src].edge.row_perm = row_perm;
    nodes_[src].edge.col_perm = col_perm;
    dirty_ = true;
}

void Graph::set_node_cost(uint32_t id, uint32_t cnot_cost) {
    if (id < nodes_.size()) {
        nodes_[id].cnot_cost = cnot_cost;
        dirty_ = true;
    }
}

Circuit Graph::extract_circuit(uint32_t start_id,
                                ParityMatrix const& input_pmat) const {
    Circuit circuit;
    circuit.set_num_qubits(n_);
    if (n_ == 0) return circuit;

    uint32_t current_id = start_id;
    std::unordered_set<uint32_t> visited;
    ParityMatrix pmat = input_pmat;

    while (current_id != UINT32_MAX && !visited.count(current_id) &&
           visited.size() < 1000) {
        visited.insert(current_id);
        GraphNode const* node = get_node(current_id);
        if (!node || node->edge.next_id == UINT32_MAX) break;

        if (node->edge.is_phase_gate) {
            if (node->edge.phase_col < pmat.m() &&
                pmat.is_column_one_hot(node->edge.phase_col) &&
                pmat.get_one_hot_row(node->edge.phase_col) == node->edge.phase_row) {
                int coeff = pmat.remove_column(node->edge.phase_col);
                circuit.append_phase_gates(node->edge.phase_row, coeff);
            }
        } else {
            if (node->edge.control != node->edge.target &&
                node->edge.control < pmat.n() && node->edge.target < pmat.n()) {
                circuit.cx(node->edge.control, node->edge.target);
                pmat.apply_cnot(node->edge.control, node->edge.target);
            }
        }
        current_id = node->edge.next_id;
    }
    return circuit;
}

void Graph::load(std::string const& path) {
    std::ifstream ifs(path);
    if (!ifs) return;

    std::string line, keyword;
    bool n_set = false;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (!n_set && line.find("n=") != std::string::npos) {
                size_t pos = line.find("n=");
                std::istringstream header(line.substr(pos + 2));
                uint32_t parsed_n;
                if (header >> parsed_n) { n_ = parsed_n; n_set = true; }
            }
            continue;
        }
        std::istringstream iss(line);
        if (!(iss >> keyword)) continue;

        if (keyword == "node") {
            uint32_t node_id;
            if (!(iss >> node_id)) continue;

            if (!std::getline(ifs, line)) return;
            std::istringstream dim_line(line);
            std::string dim_kw;
            uint32_t nn, mm;
            if (!(dim_line >> dim_kw >> nn >> mm) || dim_kw != "dim") continue;
            if (nn > 100 || mm > 100) continue;
            if (!n_set) { n_ = nn; n_set = true; }

            ParityMatrix pmat(nn, mm);
            std::vector<int> load_coeffs(mm);

            uint32_t cnot_cost = UINT32_MAX;
            if (!std::getline(ifs, line)) return;
            std::istringstream cost_line(line);
            std::string cost_kw;
            if ((cost_line >> cost_kw) && cost_kw == "cnot_cost")
                cost_line >> cnot_cost;

            if (!std::getline(ifs, line)) return;
            std::istringstream coeff_line(line);
            std::string coeff_kw;
            if (!(coeff_line >> coeff_kw) || coeff_kw != "coeffs") continue;
            for (uint32_t j = 0; j < mm; ++j)
                if (!(coeff_line >> load_coeffs[j])) return;
            for (uint32_t j = 0; j < mm; ++j)
                pmat.coeff(j) = load_coeffs[j];

            if (!std::getline(ifs, line) || line != "matrix") continue;
            for (uint32_t r = 0; r < nn; ++r) {
                if (!std::getline(ifs, line)) return;
                uint32_t col_count = std::min(static_cast<uint32_t>(line.size()), mm + nn);
                for (uint32_t c = 0; c < col_count; ++c)
                    if (line[c] == '1') pmat.row(r).set(c);
            }

            uint32_t nid = add_node(pmat);
            if (cnot_cost != UINT32_MAX) set_node_cost(nid, cnot_cost);

            if (!std::getline(ifs, line)) return;
            std::istringstream edge_line(line);
            std::string edge_kw;
            uint32_t next, is_phase, ctrl, tgt, pr, pc;
            if (!(edge_line >> edge_kw >> next >> is_phase >> ctrl >> tgt >> pr >> pc) ||
                edge_kw != "edge")
                continue;

            if (nid < nodes_.size()) {
                nodes_[nid].edge.next_id = next;
                nodes_[nid].edge.is_phase_gate = (is_phase != 0);
                nodes_[nid].edge.control = ctrl;
                nodes_[nid].edge.target = tgt;
                nodes_[nid].edge.phase_row = pr;
                nodes_[nid].edge.phase_col = pc;
                if (next != UINT32_MAX) {
                    uint32_t rps;
                    if (edge_line >> rps) {
                        nodes_[nid].edge.row_perm.resize(rps);
                        for (uint32_t i = 0; i < rps; ++i)
                            edge_line >> nodes_[nid].edge.row_perm[i];
                    }
                    uint32_t cps;
                    if (edge_line >> cps) {
                        nodes_[nid].edge.col_perm.resize(cps);
                        for (uint32_t i = 0; i < cps; ++i)
                            edge_line >> nodes_[nid].edge.col_perm[i];
                    }
                }
            }
        }
    }
    dirty_ = false;
}

void Graph::save(std::string const& path) const {
    std::ofstream ofs(path);
    if (!ofs) return;

    ofs << "# Parity Graph Cache for n=" << n_ << " (" << nodes_.size() << " nodes)\n\n";

    for (uint32_t id = 0; id < nodes_.size(); ++id) {
        GraphNode const& node = nodes_[id];
        ofs << "node " << id << "\n";
        ofs << "dim " << node.matrix->n() << " " << node.matrix->m() << "\n";
        ofs << "cnot_cost " << node.cnot_cost << "\n";
        ofs << "coeffs";
        for (uint32_t j = 0; j < node.matrix->m(); ++j)
            ofs << " " << node.matrix->coeffs()[j];
        ofs << "\n";
        ofs << "matrix\n";
        uint32_t cc = node.matrix->row(0).size();
        for (uint32_t r = 0; r < node.matrix->n(); ++r) {
            for (uint32_t c = 0; c < cc; ++c)
                ofs << (node.matrix->row(r).get(c) ? "1" : ".");
            ofs << "\n";
        }
        ofs << "edge " << node.edge.next_id << " "
            << (node.edge.is_phase_gate ? 1 : 0) << " "
            << node.edge.control << " " << node.edge.target << " "
            << node.edge.phase_row << " " << node.edge.phase_col;
        if (node.edge.next_id != UINT32_MAX) {
            ofs << " " << node.edge.row_perm.size();
            for (uint32_t r : node.edge.row_perm) ofs << " " << r;
            ofs << " " << node.edge.col_perm.size();
            for (uint32_t c : node.edge.col_perm) ofs << " " << c;
        }
        ofs << "\n\n";
    }
    dirty_ = false;
}

void Graph::load_file(std::string const& base_path) {
    for (uint32_t nn = 3; nn <= 10; ++nn) {
        std::string path = base_path + "_n" + std::to_string(nn) + ".txt";
        if (!graphs_[nn]) graphs_[nn] = std::make_unique<Graph>(nn);
        graphs_[nn]->set_n(nn);
        graphs_[nn]->load(path);
    }
}

void Graph::save_file(std::string const& base_path) {
    for (auto& [nn, graph] : graphs_) {
        if (!graph || graph->node_count() == 0) continue;
        std::string path = base_path + "_n" + std::to_string(nn) + ".txt";
        graph->save(path);
    }
}

Graph& Graph::get_graph(uint32_t nn) {
    if (graphs_.find(nn) == graphs_.end())
        graphs_[nn] = std::make_unique<Graph>(nn);
    return *graphs_[nn];
}

Graph const& Graph::get_graph(uint32_t nn) const {
    auto it = graphs_.find(nn);
    if (it == graphs_.end()) {
        static Graph empty_graph;
        return empty_graph;
    }
    return *it->second;
}

bool Graph::is_any_dirty() const {
    for (auto const& [nn, graph] : graphs_)
        if (graph && graph->is_dirty()) return true;
    return false;
}

void Graph::set_all_clean() {
    for (auto& [nn, graph] : graphs_)
        if (graph) graph->set_clean();
}

}  // namespace exact_t
