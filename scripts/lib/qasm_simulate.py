"""Load a QASM (output of exact-t-synth), simulate via qiskit-aer state vector,
return truth table mapping basis-state inputs to basis-state outputs.
"""
import re
from itertools import product
from typing import Tuple, List, Dict


def parse_qasm_io(qasm_path: str) -> Tuple[List[int], List[int]]:
    """Reads // input_qubits: ... and // output_qubits: ... header comments."""
    inputs, outputs = None, None
    with open(qasm_path) as f:
        for line in f:
            m = re.match(r'//\s*input_qubits:\s*(.+)', line)
            if m:
                inputs = [int(x) for x in m.group(1).split(',')]
            m = re.match(r'//\s*output_qubits:\s*(.+)', line)
            if m:
                outputs = [int(x) for x in m.group(1).split(',')]
            if inputs is not None and outputs is not None:
                break
    if inputs is None or outputs is None:
        raise ValueError(f"missing input_qubits/output_qubits comment in {qasm_path}")
    return inputs, outputs


def run_all_inputs(qasm_path: str,
                    input_qubits: List[int],
                    output_qubits: List[int]) -> Dict[Tuple[int, ...], Tuple[int, ...]]:
    """For each 2^|inputs| combination, prepend X gates to set basis state,
    run state-vector sim, read max-amplitude basis state, extract output bits.
    Asserts |amplitude|^2 > 0.99 on every case (deterministic Boolean output)."""
    import qiskit
    from qiskit import qasm2
    import qiskit_aer

    with open(qasm_path) as f:
        qasm_text = f.read()

    base_circuit = qasm2.loads(qasm_text)
    n_qubits = base_circuit.num_qubits
    sim = qiskit_aer.AerSimulator(method='statevector')

    table: Dict[Tuple[int, ...], Tuple[int, ...]] = {}

    for bits in product((0, 1), repeat=len(input_qubits)):
        circ = qiskit.QuantumCircuit(n_qubits)
        for q, b in zip(input_qubits, bits):
            if b == 1:
                circ.x(q)
        circ.compose(base_circuit, inplace=True)
        circ.save_statevector()

        result = sim.run(qiskit.transpile(circ, sim)).result()
        sv = result.get_statevector(circ)

        amps = abs(sv.data) ** 2
        argmax_idx = int(amps.argmax())
        max_prob = float(amps[argmax_idx])
        if max_prob < 0.99:
            raise RuntimeError(
                f"non-deterministic output at input {bits}: "
                f"max basis prob = {max_prob:.4f} (expected ~1.0)"
            )

        out_bits = tuple(
            (argmax_idx >> q) & 1 for q in output_qubits
        )
        table[bits] = out_bits

    return table
