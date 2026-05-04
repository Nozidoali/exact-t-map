#include "core/circuit.hpp"

#include <algorithm>

namespace exact_t {

uint32_t Circuit::num_qubits() const { return num_qubits_; }
void Circuit::set_num_qubits(uint32_t n) { num_qubits_ = n; }

std::vector<Gate> const& Circuit::gates() const { return gates_; }

std::vector<uint32_t> const& Circuit::inputs() const { return inputs_; }
std::vector<uint32_t>& Circuit::inputs() { return inputs_; }

std::vector<uint32_t> const& Circuit::outputs() const { return outputs_; }
std::vector<uint32_t>& Circuit::outputs() { return outputs_; }

void Circuit::add_input(uint32_t q) { inputs_.push_back(q); }
void Circuit::add_output(uint32_t q) { outputs_.push_back(q); }

void Circuit::h(uint32_t q) {
    gates_.push_back({GateType::H, q});
}

void Circuit::s(uint32_t q) {
    gates_.push_back({GateType::S, q});
}

void Circuit::sdg(uint32_t q) {
    gates_.push_back({GateType::Sdg, q});
}

void Circuit::t(uint32_t q) {
    gates_.push_back({GateType::T, q});
}

void Circuit::tdg(uint32_t q) {
    gates_.push_back({GateType::Tdg, q});
}

void Circuit::cx(uint32_t ctrl, uint32_t tgt) {
    gates_.push_back({GateType::CX, tgt, ctrl});
}

void Circuit::z(uint32_t q) {
    gates_.push_back({GateType::Z, q});
}

void Circuit::x(uint32_t q) {
    gates_.push_back({GateType::X, q});
}

void Circuit::append_phase_gates(uint32_t qubit, int coeff) {
    coeff = ((coeff % 8) + 8) % 8;
    switch (coeff) {
        case 0: break;
        case 1: t(qubit); break;
        case 2: s(qubit); break;
        case 3: s(qubit); t(qubit); break;
        case 4: z(qubit); break;
        case 5: z(qubit); t(qubit); break;
        case 6: sdg(qubit); break;
        case 7: tdg(qubit); break;
    }
}

void Circuit::append(Circuit const& other, std::vector<uint32_t> const& mapping) {
    for (auto const& gate : other.gates_) {
        Gate mapped_gate = gate;
        mapped_gate.target = mapping[gate.target];
        if (gate.type == GateType::CX)
            mapped_gate.control = mapping[gate.control];
        gates_.push_back(mapped_gate);
    }
}

void Circuit::toffoli(uint32_t ctrl1, uint32_t ctrl2, uint32_t target) {
    h(target);
    cx(ctrl2, target); tdg(target);
    cx(ctrl1, target); t(target);
    cx(ctrl2, target); tdg(target);
    cx(ctrl1, target);
    t(target); h(target);
    cx(ctrl1, ctrl2); tdg(ctrl2);
    cx(ctrl1, ctrl2);
    t(ctrl1); t(ctrl2);
}

void Circuit::finalize_qubit_count() {
    uint32_t max_q = 0;
    for (auto const& gate : gates_) {
        max_q = std::max(max_q, gate.target);
        if (gate.type == GateType::CX)
            max_q = std::max(max_q, gate.control);
    }
    for (uint32_t q : inputs_) max_q = std::max(max_q, q);
    for (uint32_t q : outputs_) max_q = std::max(max_q, q);
    num_qubits_ = max_q + 1;
}

GateCounts Circuit::count_gates() const {
    GateCounts gc;
    for (auto const& g : gates_) {
        switch (g.type) {
            case GateType::T: case GateType::Tdg: ++gc.t; break;
            case GateType::CX: ++gc.cx; break;
            case GateType::H: ++gc.h; break;
            case GateType::S: case GateType::Sdg: ++gc.s; break;
            case GateType::Z: ++gc.z; break;
            case GateType::X: ++gc.x; break;
        }
    }
    return gc;
}

}  // namespace exact_t
