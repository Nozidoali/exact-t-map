/*! \file solver.hpp
 *  \brief CliffordTSolver: exact T-count synthesis from tensor targets.
 *
 *  Finds minimum T-count phase polynomials for Boolean functions
 *  using Gaussian elimination over GF(2) with basis refinement.
 *
 *  Algorithm:
 *  1. Lazily create and cache SolverTables per-n on first use
 *  2. Greedy CNOT preprocessing to reduce tensor complexity
 *  3. Gaussian elimination for initial phase polynomial
 *  4. Coordinate descent basis refinement to minimize T-count
 *  5. Z8 coefficient correction for exact circuit synthesis
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/bits.hpp"
#include "core/circuit.hpp"
#include "core/monomial.hpp"
#include "core/phase.hpp"
#include "core/tensor.hpp"

namespace exact_t {

class BilinearFunction;
class SolverTables;

/*! \brief Phase polynomial with mod-8 coefficients.
 *
 *  Each entry maps a parity mask to a Z8 coefficient representing
 *  the number of T-gate rotations (mod 8) on that parity.
 */
struct Z8Phase {
    uint32_t n{0};
    std::unordered_map<uint32_t, int> terms;

    Z8Phase();
    explicit Z8Phase(uint32_t n);

    void set_coeff(uint32_t parity, int coeff);
    int get_coeff(uint32_t parity) const;
    void add_coeff(uint32_t parity, int val);
    bool empty() const;
    uint32_t size() const;
    std::vector<std::pair<uint32_t, int>> to_vector() const;

    /*! \brief Convert mod-2 Phase to Z8Phase (each term gets coeff 1). */
    static Z8Phase from_phase(Phase const& phase);

    Z8Phase operator+(Z8Phase const& o) const;

    /*! \brief Filter terms by qubit involvement.
     *  \return Pair of (terms involving output qubits or global, terms not involving them)
     */
    std::pair<Z8Phase, Z8Phase> filter_by_qubits(uint32_t qubit_mask) const;
};

/*! \brief Precomputed Gaussian elimination transform. */
struct Transform {
    std::vector<Bits> matrix;
    std::vector<int> row_to_cand;
    uint32_t rank{0};
};

/*! \brief Solver configuration. */
struct SolverConfig {
    uint32_t max_flip_depth{4};  /*!< Max simultaneous basis flips (1-4) */
};

/*! \brief Clifford+T exact synthesis solver.
 *
 *  Accepts a maximum n at construction. SolverTables are lazily created
 *  and cached per-n on first use, so small cuts avoid the cost of large
 *  table precomputation. synthesize() auto-detects n from the input.
 */
class CliffordTSolver {
  public:
    CliffordTSolver(uint32_t n, SolverConfig const& config = {});

    /*! \brief Returns max_n (upper bound passed at construction). */
    uint32_t n() const { return max_n_; }
    SolverConfig const& config() const { return config_; }

    /*! \brief Find optimal phase polynomial for a tensor target.
     *
     *  If the tensor has don't-care bits set, they provide extra
     *  freedom during basis refinement, potentially reducing T-count.
     */
    Phase solve(Tensor const& target) const;

    /*! \brief Full pipeline: BilinearFunction -> optimized Clifford+T circuit.
     *
     *  Auto-detects n from bf.total_qubits().
     *  1. Convert to tensor and solve for optimal phase polynomial
     *  2. Apply Z8 coefficient correction
     *  3. Synthesize circuit with CNOT chains and T/S/Z rotations
     */
    Circuit synthesize(BilinearFunction const& bf) const;

    /*! \brief Synthesize with optional clean-ancilla don't-care optimization.
     *
     *  When use_ancilla_dc is true, uses max_n_ to preserve ancilla
     *  padding-as-don't-care optimization via a two-level approach:
     *  1. Kernel augmentation: augments basis with output kernel vectors
     *  2. Z8 coefficient optimization: tries all compensation combinations
     *     to minimize T-count via don't-care monomials
     */
    Circuit synthesize(BilinearFunction const& bf, bool use_ancilla_dc) const;

    /*! \brief Synthesize from tensor with optional don't-care mask.
     *
     *  Uses n from target.num_vars(). When don't-cares are present,
     *  the solver exploits them to find a lower T-count phase polynomial.
     */
    Circuit synthesize_tensor(Tensor const& target,
                               uint32_t num_inputs, uint32_t num_outputs) const;

    /*! \brief Load solver cache from CSV file (hlqcs format). */
    void load_cache(std::string const& path);

    /*! \brief Number of cached solutions. */
    size_t cache_size() const { return cache_.size(); }

    /*! \brief Number of cache hits since construction. */
    size_t cache_hits() const { return cache_hits_; }

    /*! \brief Number of cache misses since construction. */
    size_t cache_misses() const { return cache_misses_; }

private:
    uint32_t max_n_;
    SolverConfig config_;
    mutable std::unordered_map<uint32_t, std::shared_ptr<SolverTables>> tables_cache_;

    mutable std::unordered_map<Tensor, Phase> cache_;
    mutable size_t cache_hits_{0};
    mutable size_t cache_misses_{0};

    /*! \brief Get or lazily create SolverTables for the given n. */
    std::shared_ptr<SolverTables> get_tables(uint32_t n) const;

    Z8Monomial phase_to_z8(Phase const& phase, uint32_t n) const;
    Z8Phase solve_z8(Z8Monomial const& diff, uint32_t n) const;

    std::vector<std::pair<uint32_t, uint32_t>> greedy_cnot_preprocess(
        Tensor& target, uint32_t n) const;
    void emit_phase_rotations(Circuit& circ, Z8Phase const& z8_phase,
                               uint32_t n) const;

    /*! \brief Solve with additional kernel vectors for basis augmentation. */
    Phase solve_with_kernels(Tensor const& target,
                              std::vector<Bits> const& extra_kernels) const;
};

/*! \brief Full Clifford+T synthesis: BilinearFunction -> optimized circuit.
 *
 *  Convenience free function wrapping CliffordTSolver::synthesize.
 */
Circuit synthesize_clifford_t(CliffordTSolver const& solver,
                               BilinearFunction const& bf);

}  // namespace exact_t
