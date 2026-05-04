# Quantum Circuit Synthesis Using an Exact T Library

[![CI](https://github.com/Nozidoali/exact-t-map/actions/workflows/ci.yml/badge.svg)](https://github.com/Nozidoali/exact-t-map/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![CMake](https://img.shields.io/badge/CMake-%E2%89%A53.16-064F8C?logo=cmake&logoColor=white)](https://cmake.org/)
[![DAC '26](https://img.shields.io/badge/DAC-2026-red.svg)](#citing)

Exact-synthesis-driven Clifford+T technology mapper for XAG (XOR-AND graph)
networks. Lowers a Boolean network into a Clifford+T circuit with provably
minimum per-cut T-count via bilinear cut enumeration and an exact
phase-polynomial solver.

## Build

```bash
cmake -B build && cmake --build build -j
```

Produces `build/exact-t-synth`. Requires C++17 and CMake $\ge$ 3.16.

## Usage

```bash
./build/exact-t-synth -n 8 --cut-size 6 \
    -o multiplier.qasm benchmarks/epfl/multiplier.v
```

Common flags: `--cut-size`, `--max-cuts-per-node`, `--stats-out`. See
`./build/exact-t-synth --help` for the full list. Try a small test first:

```bash
./build/exact-t-synth -n 5 --cut-size 4 \
    -o maj3.qasm benchmarks/correctness/maj3.v
```

## Verify

```bash
python3 scripts/verify-correctness.py benchmarks/correctness/*.v
```

Requires `qiskit` and `qiskit-aer`. Each test maps a Verilog file to QASM,
simulates the circuit end-to-end, and compares against the truth table
evaluated directly from the Verilog.

## Repository

```
src/{core,network,mapper,solver,io}    C++ implementation
scripts/                                Correctness verification harness
benchmarks/correctness/                 Small Verilog test cases
```

## Citing

```bibtex
@inproceedings{WangDAC2026ExactT,
  author    = {Wang, Hanyu and Yu, Mingfei and Wu, Xinrui and Cong, Jason},
  title     = {Quantum Circuit Synthesis Using an Exact T Library},
  booktitle = {Proceedings of the 63rd ACM/IEEE Design Automation Conference (DAC '26)},
  year      = {2026},
  address   = {Long Beach, CA, USA},
  month     = jul,
  note      = {July 26--29, 2026}
}
```

## License

MIT. See [LICENSE](LICENSE).
