/*! \file tensor.hpp
 *  \brief Tensor: binary monomial representation over n variables.
 *
 *  Each bit position (0 to 2^n-1) corresponds to a monomial mask.
 *  Bit at position m is set iff the monomial prod_{i in m} x_i is present.
 *  Used as the solver interface for Clifford+T synthesis.
 */

#pragma once

#include <cstdint>
#include <functional>
#include <string>

#include "core/bits.hpp"

namespace exact_t {

class Phase;
class Z8Monomial;

/*! \brief Binary monomial representation with 2^n entries.
 *
 *  Optional care mask: when set, only care-bit positions constrain
 *  the solver. Don't-care positions give extra synthesis freedom.
 */
class Tensor {
public:
    Tensor();
    explicit Tensor(uint32_t n);

    /*! \brief Number of variables. */
    uint32_t num_vars() const;

    /*! \brief Read-only access to data bits. */
    Bits const& data() const;

    /*! \brief Mutable access to data bits. */
    Bits& data();

    /*! \brief Read-only access to care mask. */
    Bits const& care() const;

    /*! \brief Mutable access to care mask. */
    Bits& care();

    /*! \brief Toggle linear monomial x_i. */
    void add_linear(uint32_t i);

    /*! \brief Toggle quadratic monomial x_i * x_j. */
    void add_quadratic(uint32_t i, uint32_t j);

    /*! \brief Toggle cubic monomial x_i * x_j * x_k. */
    void add_cubic(uint32_t i, uint32_t j, uint32_t k);

    /*! \brief Basis change: x_target -> x_target XOR x_ctrl.
     *
     *  For each monomial with target bit set, toggles the ctrl bit.
     */
    void apply_cnot(uint32_t ctrl, uint32_t target);

    /*! \brief Mark monomial at position mask as don't-care. */
    void set_dont_care(uint32_t mask);

    /*! \brief Mark all monomials involving variable i as don't-care. */
    void set_output_dont_care(uint32_t output_bit);

    /*! \brief Check if position has a care constraint. */
    bool is_care(uint32_t mask) const;

    /*! \brief True if any don't-care bits are set. */
    bool has_dont_cares() const;

    /*! \brief Number of don't-care positions. */
    uint32_t num_dont_cares() const;

    /*! \brief Count of nonzero monomials. */
    uint32_t num_terms() const;

    /*! \brief Binary string of monomial presence bits. */
    std::string to_string() const;

    Tensor& operator+=(Tensor const& o);
    Tensor operator+(Tensor const& o) const;
    bool operator==(Tensor const& o) const;

    /*! \brief Build tensor from all 3-wise products of bits in mask.
     *  \param n Number of variables
     *  \param mask Bitmask selecting variables
     */
    static Tensor from_parity_mask(uint32_t n, uint32_t mask);

    /*! \brief Convert to Z8 monomial representation.
     *
     *  Each set monomial bit gets coefficient = t_coefficient(degree).
     */
    Z8Monomial to_z8_monomial() const;

    /*! \brief Build tensor from phase polynomial (XOR of per-term tensors). */
    static Tensor from_phase(Phase const& p);

private:
    uint32_t n_{0};
    Bits data_;
    Bits care_;  /*!< Care mask (empty = all care). */
};

}  // namespace exact_t

namespace std {
template <>
struct hash<exact_t::Tensor> {
    size_t operator()(exact_t::Tensor const& t) const {
        size_t h = 14695981039346656037ULL;
        for (uint64_t w : t.data().raw()) {
            h ^= w;
            h *= 1099511628211ULL;
        }
        if (t.care().size() > 0) {
            for (uint64_t w : t.care().raw()) {
                h ^= w;
                h *= 1099511628211ULL;
            }
        }
        return h;
    }
};
}  // namespace std
