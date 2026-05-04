cmake_minimum_required(VERSION 3.16)
project(exact-t LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Core library (no external dependencies)
add_library(exact-t-core
    src/core/bits.cpp
    src/core/circuit.cpp
    src/core/tensor.cpp
    src/core/phase.cpp
    src/core/monomial.cpp
    src/core/bilinear.cpp
    src/solver/solver.cpp
    src/solver/tables.cpp
    src/solver/gauss.cpp
    src/solver/closed_form.cpp
    src/io/circuit-io.cpp
    src/io/trace.cpp
    src/io/stats.cpp
    src/solver/parity-matrix.cpp
    src/solver/graph.cpp
    src/solver/parity-solver.cpp
    src/solver/affine-state.cpp
    src/solver/affine.cpp
)
target_include_directories(exact-t-core PUBLIC src)

# Full library (with XAG network and mapper)
add_library(exact-t
    src/network/xag.cpp
    src/network/reader.cpp
    src/mapper/cut.cpp
    src/mapper/allocate.cpp
    src/mapper/bilinear.cpp
    src/mapper/evaluator.cpp
    src/mapper/daomap.cpp
    src/mapper/extractor.cpp
    src/mapper/mapper.cpp
    src/io/io.cpp
    src/io/extraction-dump.cpp
)
target_include_directories(exact-t PUBLIC src)
target_link_libraries(exact-t PUBLIC exact-t-core)

# CLI binary
add_executable(exact-t-synth src/main.cpp)
target_link_libraries(exact-t-synth exact-t)
