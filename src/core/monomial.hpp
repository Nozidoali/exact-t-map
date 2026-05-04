/*! \file monomial.hpp
 *  \brief Z8Monomial: polynomial with mod-8 coefficients indexed by monomial masks.
 *
 *  Each coefficient at index mask represents the monomial prod_{i in mask} x_i.
 *  Used in the Clifford+T solver for exact phase arithmetic.
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace exact_t {

/*! \brief Polynomial with Z8 coefficients indexed by monomial masks.
 *
 *  Stores coefficients mod 8 for each possible monomial up to degree n.
 *  The coefficient at index mask represents the monomial prod_{i in mask} x_i.
 */
class Z8Monomial {
public:
    Z8Monomial();
    explicit Z8Monomial(uint32_t n);

    /*! \brief Number of variables. */
    uint32_t num_vars() const;

    /*! \brief Read-only access to coefficients. */
    std::vector<int> const& coeffs() const;

    /*! \brief Normalize value to unsigned mod-8 range [0, 7]. */
    static int mod8(int val) { return ((val % 8) + 8) % 8; }

    /*! \brief Convert unsigned mod-8 value to signed range [-4, 3]. */
    static int to_signed(int val) { return (val > 4) ? val - 8 : val; }

    /*! \brief Monomial degree from mask (popcount). */
    static uint32_t degree(uint32_t mask) {
        return static_cast<uint32_t>(__builtin_popcount(mask));
    }

    /*! \brief T-gate coefficient for monomial degree.
     *
     *  Based on inclusion-exclusion for exp(i*pi/4 * prod(z_i)):
     *  - Degree 0-1: coefficient 1
     *  - Degree 2: coefficient 6 (= -2 mod 8)
     *  - Degree 3+: coefficient 4
     */
    static int t_coefficient(uint32_t mask) {
        uint32_t deg = degree(mask);
        return (deg <= 1) ? 1 : (deg == 2) ? 6 : 4;
    }

    /*! \brief Divisor for monomial degree normalization.
     *
     *  Returns 2^(deg-1) for deg > 1, else 1. Used in
     *  solve_linear_system_z8 to invert the parity-to-monomial map.
     */
    static int degree_divisor(uint32_t mask) {
        uint32_t deg = degree(mask);
        return (deg <= 1) ? 1 : (1 << (deg - 1));
    }

    /*! \brief Check if coefficient is divisible by 2*t_coefficient. */
    static bool is_divisible(int coeff, uint32_t mask) {
        int signed_coeff = to_signed(mod8(coeff));
        return (signed_coeff % (t_coefficient(mask) * 2)) == 0;
    }

    void add_coeff(uint32_t idx, int val);

    int& operator[](uint32_t idx) { return coeffs_[idx]; }
    int const& operator[](uint32_t idx) const { return coeffs_[idx]; }
    size_t size() const { return coeffs_.size(); }

    /*! \brief Check if all coefficients are even. */
    bool is_even() const;

    Z8Monomial operator+(Z8Monomial const& o) const;
    Z8Monomial operator-(Z8Monomial const& o) const;
    Z8Monomial operator*(int scalar) const;

    std::string to_string() const;

private:
    uint32_t n_{0};
    std::vector<int> coeffs_;
};

}  // namespace exact_t
