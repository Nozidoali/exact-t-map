/*! \file bilinear.hpp
 *  \brief BilinearFunction: Boolean functions as sums of bilinear terms.
 *
 *  Encodes f_k = XOR_{(i,j)} (x_i AND x_j). Each term is a triplet
 *  (i, j, k) meaning "input i AND input j contributes to output k".
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "core/bits.hpp"
#include "core/circuit.hpp"
#include "core/monomial.hpp"
#include "core/phase.hpp"
#include "core/tensor.hpp"

namespace exact_t {

/*! \brief A single bilinear term: x_i AND x_j -> output k. */
struct Triplet {
    uint8_t i;  /*!< First input index */
    uint8_t j;  /*!< Second input index */
    uint8_t k;  /*!< Output index */
};

/*! \brief Boolean function represented as a sum of bilinear terms.
 *
 *  Each output bit is the XOR of AND terms over pairs of input bits.
 *  Provides two circuit synthesis paths:
 *  - to_circuit(): one Toffoli per triplet (7T each, simple)
 *  - to_phase_circuit(): phase polynomial synthesis with term cancellation
 */
class BilinearFunction {
public:
    BilinearFunction(uint32_t ni, uint32_t no);

    /*! \brief Number of input variables. */
    uint32_t num_inputs() const;

    /*! \brief Number of output variables. */
    uint32_t num_outputs() const;

    /*! \brief Read-only access to triplets. */
    std::vector<Triplet> const& triplets() const;

    /*! \brief Add a bilinear term: x_i AND x_j contributes to output k. */
    void add_triplet(uint32_t i, uint32_t j, uint32_t k);

    /*! \brief Total qubits needed: num_inputs + num_outputs. */
    uint32_t total_qubits() const;

    /*! \brief Human-readable representation. */
    std::string to_string() const;

    /*! \brief Convert to cubic monomial tensor representation.
     *
     *  Each triplet (i, j, k) becomes cubic monomial
     *  x_i * x_j * x_{num_inputs+k} in a tensor of dimension dim.
     */
    Tensor to_tensor(uint32_t dim) const;

    /*! \brief Convert to tensor with dim = total_qubits(). */
    Tensor to_tensor() const;

    /*! \brief Collect parity bitmasks for triplets targeting a given output. */
    std::vector<uint32_t> get_output_parities(uint32_t output_idx) const;

    /*! \brief Get kernel vectors for each output (OR-closure of parities).
     *
     *  For each output, computes the BFS OR-closure of its input parities,
     *  then XORs the corresponding parity-to-monomial vectors to produce
     *  a kernel vector for basis augmentation.
     */
    std::vector<Bits> get_kernels(uint32_t dim) const;

    /*! \brief Get don't-care monomials for Z8 coefficient optimization.
     *
     *  For each output, computes the inclusion-exclusion polynomial
     *  over its input parities, yielding a Z8Monomial that can absorb
     *  arbitrary phase rotations via compensation gates.
     */
    std::vector<Z8Monomial> get_dontcare_monomials(uint32_t dim) const;

    /*! \brief Convert to phase polynomial via CCZ expansion.
     *
     *  Each triplet generates 7 mod-2 phase terms from the
     *  CCZ inclusion-exclusion formula. Terms from different
     *  triplets may cancel.
     */
    Phase to_phase() const;

    /*! \brief Decompose into Clifford+T via individual Toffoli gates.
     *
     *  Allocates num_inputs + num_outputs qubits. For each triplet
     *  (i, j, k), applies a Toffoli on (i, j, num_inputs + k)
     *  using the standard 7-T-gate decomposition.
     */
    Circuit to_circuit() const;

    /*! \brief Decompose into Clifford+T via phase polynomial synthesis.
     *
     *  Expands all triplets as CCZ gates into a single phase polynomial,
     *  allowing term cancellation to reduce T-count. Synthesizes each
     *  remaining term as CNOT-parity + rotation + uncompute.
     *
     *  Uses mod-8 coefficient tracking: T^2=S, T^4=Z, T^6=Sdg, T^7=Tdg
     *  so non-T rotations become Clifford gates.
     */
    Circuit to_phase_circuit() const;

private:
    uint32_t num_inputs_;
    uint32_t num_outputs_;
    std::vector<Triplet> triplets_;
};

/*! \brief Convert to tensor with ancilla qubits marked as don't-care.
 *
 *  Calls bf.to_tensor(dim), then marks qubits [bf.total_qubits(), dim)
 *  as output don't-cares, giving the solver extra freedom.
 */
Tensor to_tensor_with_ancilla_dc(BilinearFunction const& bf, uint32_t dim);

}  // namespace exact_t
