#!/usr/bin/env python3
"""Verify exact-t-synth's QASM output is functionally equivalent to input Verilog.

Usage:
    scripts/verify-correctness.py [--verbose] [--keep-qasm] benchmarks/correctness/*.v

Each .v is mapped via the exact-t-synth binary (under SLURM), then the resulting
QASM is simulated for all 2^|inputs| input combinations and compared against the
truth table directly evaluated from the Verilog.

Exit codes:
    0  all tests pass
    1  one or more tests fail
    2  qiskit not installed
"""
import argparse
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'scripts'))

from lib import verilog_truth  # noqa: E402


def map_verilog(v_path: str, qasm_path: str, verbose: bool):
    cmd = [
        'srun', '-p', 'batch', '--qos=low', '--time=00:01:00',
        os.path.join(REPO, 'build', 'exact-t-synth'),
        '-o', qasm_path, v_path,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(f"mapper FAILED on {v_path} (exit {proc.returncode})\n")
        sys.stderr.write(proc.stderr)
        return False
    if verbose and proc.stderr:
        sys.stderr.write(f"--- mapper stderr for {v_path} ---\n{proc.stderr}")
    return True


def verify_one(v_path: str, verbose: bool, keep_qasm: bool) -> bool:
    name = os.path.splitext(os.path.basename(v_path))[0]
    qasm_path = os.path.join(REPO, 'results', 'correctness', f'{name}.qasm')
    os.makedirs(os.path.dirname(qasm_path), exist_ok=True)

    if not map_verilog(v_path, qasm_path, verbose):
        print(f"{name:<12s}  FAIL  mapper failure")
        return False

    expected = verilog_truth.compute_truth_table(v_path)
    in_names, out_names = verilog_truth.io_signals(v_path)

    from lib import qasm_simulate
    in_q, out_q = qasm_simulate.parse_qasm_io(qasm_path)

    if len(in_q) != len(in_names):
        print(f"{name:<12s}  FAIL  qasm has {len(in_q)} input qubits, "
              f"verilog has {len(in_names)} inputs")
        return False
    if len(out_q) != len(out_names):
        print(f"{name:<12s}  FAIL  qasm has {len(out_q)} output qubits, "
              f"verilog has {len(out_names)} outputs")
        return False

    actual = qasm_simulate.run_all_inputs(qasm_path, in_q, out_q)

    cases = len(expected)
    for inp, exp_out in sorted(expected.items()):
        act_out = actual.get(inp)
        if act_out != exp_out:
            print(f"{name:<12s}  FAIL  first mismatch at input={inp}: "
                  f"expected={exp_out} actual={act_out}")
            if verbose:
                bits_diff = [i for i, (e, a) in enumerate(zip(exp_out, act_out)) if e != a]
                print(f"               differing output bits: {bits_diff}")
            return False

    print(f"{name:<12s}  PASS  ({len(in_names)} inputs, {cases} cases)")

    if not keep_qasm:
        try:
            os.remove(qasm_path)
        except OSError:
            pass

    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--keep-qasm', action='store_true')
    parser.add_argument('verilog', nargs='+', help='Verilog files to verify')
    args = parser.parse_args()

    try:
        import qiskit  # noqa: F401
        import qiskit_aer  # noqa: F401
    except ImportError:
        print("qiskit not installed. Run:", file=sys.stderr)
        print("  python3 -m pip install --user qiskit qiskit-aer", file=sys.stderr)
        return 2

    all_ok = True
    for v in args.verilog:
        ok = verify_one(v, args.verbose, args.keep_qasm)
        all_ok = all_ok and ok

    if all_ok:
        print(f"\nAll {len(args.verilog)} tests passed.")
        return 0
    else:
        print(f"\nSome tests failed.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
