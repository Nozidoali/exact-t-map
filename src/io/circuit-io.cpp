#include "io/circuit-io.hpp"

#include <fstream>
#include <sstream>

namespace exact_t {

namespace {

char const* gate_name_qasm(GateType type) {
    switch (type) {
        case GateType::H: return "h";
        case GateType::S: return "s";
        case GateType::Sdg: return "sdg";
        case GateType::T: return "t";
        case GateType::Tdg: return "tdg";
        case GateType::Z: return "z";
        case GateType::X: return "x";
        case GateType::CX: return "cx";
    }
    return "?";
}

}  // namespace

std::string to_qasm(Circuit const& circuit) {
    std::ostringstream os;
    os << "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n";

    if (!circuit.inputs().empty()) {
        os << "// n_inputs: " << circuit.inputs().size() << "\n";
        os << "// input_qubits:";
        for (size_t i = 0; i < circuit.inputs().size(); ++i)
            os << (i ? "," : " ") << circuit.inputs()[i];
        os << "\n";
    }
    if (!circuit.outputs().empty()) {
        os << "// n_outputs: " << circuit.outputs().size() << "\n";
        os << "// output_qubits:";
        for (size_t i = 0; i < circuit.outputs().size(); ++i)
            os << (i ? "," : " ") << circuit.outputs()[i];
        os << "\n";
    }

    os << "qreg q[" << circuit.num_qubits() << "];\n";
    if (!circuit.outputs().empty())
        os << "creg c[" << circuit.outputs().size() << "];\n";
    os << "\n";

    for (auto const& gate : circuit.gates()) {
        if (gate.type == GateType::CX) {
            os << "cx q[" << gate.control << "],q[" << gate.target << "];\n";
        } else {
            os << gate_name_qasm(gate.type) << " q[" << gate.target << "];\n";
        }
    }

    return os.str();
}

std::string to_qc(Circuit const& circuit) {
    std::ostringstream os;
    os << ".v";
    for (uint32_t i = 0; i < circuit.num_qubits(); ++i)
        os << " " << i;
    os << "\n";

    if (!circuit.inputs().empty()) {
        os << ".i";
        for (uint32_t q : circuit.inputs()) os << " " << q;
        os << "\n";
    }
    if (!circuit.outputs().empty()) {
        os << ".o";
        for (uint32_t q : circuit.outputs()) os << " " << q;
        os << "\n";
    }

    for (auto const& gate : circuit.gates()) {
        switch (gate.type) {
            case GateType::H: os << "H " << gate.target << "\n"; break;
            case GateType::S: os << "S " << gate.target << "\n"; break;
            case GateType::Sdg: os << "S* " << gate.target << "\n"; break;
            case GateType::T: os << "T " << gate.target << "\n"; break;
            case GateType::Tdg: os << "T* " << gate.target << "\n"; break;
            case GateType::Z: os << "Z " << gate.target << "\n"; break;
            case GateType::X: os << "X " << gate.target << "\n"; break;
            case GateType::CX:
                os << "CNOT " << gate.control << " " << gate.target << "\n";
                break;
        }
    }
    return os.str();
}

void write_qasm(Circuit const& circuit, std::string const& path) {
    std::ofstream f(path);
    if (f.is_open()) f << to_qasm(circuit);
}

void write_qc(Circuit const& circuit, std::string const& path) {
    std::ofstream f(path);
    if (f.is_open()) f << to_qc(circuit);
}

}  // namespace exact_t
