/*! \file circuit.hpp
 *  \brief Lightweight Clifford+T circuit representation.
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace exact_t {

/*! \brief Gate type in a Clifford+T circuit. */
enum class GateType : uint8_t { H, S, Sdg, T, Tdg, CX, Z, X };

/*! \brief A single gate in a circuit.
 *
 *  For single-qubit gates, only target is used.
 *  For CX, control and target are both used.
 */
struct Gate {
    GateType type;
    uint32_t target;
    uint32_t control{0};
};

/*! \brief Gate count summary. */
struct GateCounts {
    uint32_t t{0};
    uint32_t cx{0};
    uint32_t h{0};
    uint32_t s{0};
    uint32_t z{0};
    uint32_t x{0};
};

/*! \brief A flat gate-list circuit for Clifford+T gates. */
class Circuit {
public:
    /*! \brief Number of qubits in the circuit. */
    uint32_t num_qubits() const;

    /*! \brief Set the number of qubits. */
    void set_num_qubits(uint32_t n);

    /*! \brief Read-only access to gates. */
    std::vector<Gate> const& gates() const;

    /*! \brief Read-only access to input qubits. */
    std::vector<uint32_t> const& inputs() const;

    /*! \brief Mutable access to input qubits. */
    std::vector<uint32_t>& inputs();

    /*! \brief Read-only access to output qubits. */
    std::vector<uint32_t> const& outputs() const;

    /*! \brief Mutable access to output qubits. */
    std::vector<uint32_t>& outputs();

    /*! \brief Add an input qubit. */
    void add_input(uint32_t q);

    /*! \brief Add an output qubit. */
    void add_output(uint32_t q);

    /*! \brief Append a Hadamard gate. */
    void h(uint32_t q);

    /*! \brief Append an S gate (phase pi/2). */
    void s(uint32_t q);

    /*! \brief Append an S-dagger gate (phase -pi/2). */
    void sdg(uint32_t q);

    /*! \brief Append a T gate (phase pi/4). */
    void t(uint32_t q);

    /*! \brief Append a T-dagger gate (phase -pi/4). */
    void tdg(uint32_t q);

    /*! \brief Append a CNOT gate. */
    void cx(uint32_t ctrl, uint32_t tgt);

    /*! \brief Append a Pauli-Z gate. */
    void z(uint32_t q);

    /*! \brief Append a Pauli-X gate. */
    void x(uint32_t q);

    /*! \brief Append Z-rotation gates for a Z8 phase coefficient.
     *
     *  Maps coefficient mod 8: 1=T, 2=S, 3=S+T, 4=Z, 5=Z+T, 6=Sdg, 7=Tdg.
     */
    void append_phase_gates(uint32_t qubit, int coeff);

    /*! \brief Append another circuit with qubit remapping. */
    void append(Circuit const& other, std::vector<uint32_t> const& mapping);

    /*! \brief Append Toffoli gate (7-T decomposition). */
    void toffoli(uint32_t ctrl1, uint32_t ctrl2, uint32_t target);

    /*! \brief Update num_qubits to cover all referenced qubits. */
    void finalize_qubit_count();

    /*! \brief Count gates by type. */
    GateCounts count_gates() const;

private:
    uint32_t num_qubits_{0};
    std::vector<Gate> gates_;
    std::vector<uint32_t> inputs_;
    std::vector<uint32_t> outputs_;
};

}  // namespace exact_t
