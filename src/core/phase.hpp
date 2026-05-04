/*! \file phase.hpp
 *  \brief Phase polynomial: sum of pi/4 rotations on parity terms.
 *
 *  Represents exp(i pi/4 * sum prod_{j in S} z_j) with mod-2
 *  arithmetic on parity masks. Each term is a bitmask selecting
 *  which qubits participate in the parity.
 */

#pragma once

#include <cstdint>
#include <set>
#include <string>

namespace exact_t {

/*! \brief Phase polynomial with mod-2 parity terms.
 *
 *  Each entry in the set is a bitmask: bit i set means variable z_i
 *  participates. Toggle semantics: adding a term twice cancels it.
 */
class Phase {
public:
    Phase();
    explicit Phase(uint32_t n);

    /*! \brief Number of variables. */
    uint32_t num_vars() const;

    /*! \brief Read-only access to the set of terms. */
    std::set<uint32_t> const& terms() const;

    /*! \brief Toggle term with mod-2 cancellation.
     *
     *  Adds term if absent, removes if present (XOR semantics).
     *  Masks the term to n bits.
     */
    void toggle_term(uint32_t t);

    /*! \brief Add CCZ phase contribution on qubits i, j, k.
     *
     *  Toggles the 7 inclusion-exclusion terms:
     *  (ijk), (ij), (ik), (jk), (i), (j), (k).
     */
    void ccz(uint32_t i, uint32_t j, uint32_t k);

    /*! \brief Number of active terms. */
    uint32_t size() const;

    /*! \brief Check if phase polynomial is empty. */
    bool empty() const;

    /*! \brief Check if a specific term is present. */
    bool contains(uint32_t t) const;

    /*! \brief Inverse CNOT basis transformation.
     *
     *  For each parity with target bit set, toggles the ctrl bit.
     */
    void apply_inverse_cnot(uint32_t ctrl, uint32_t target);

    Phase& operator+=(Phase const& o);

    /*! \brief Permutation-invariant hash for cache lookup.
     *
     *  Computes a hash invariant under qubit permutations by
     *  canonicalizing column order of the term-by-variable matrix.
     */
    uint64_t hash() const;

    /*! \brief Format as "z0*z1 + z2 + ...". */
    std::string to_string() const;

private:
    uint32_t n_{0};
    std::set<uint32_t> terms_;
};

}  // namespace exact_t
