/*! \file bits.hpp
 *  \brief Fixed-size bit vector with efficient bitwise operations.
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace exact_t {

/*! \brief Fixed-size bit vector stored as 64-bit words.
 *
 *  Provides efficient set/get/flip and bulk XOR operations.
 *  Used by Tensor to store 2^n monomial presence bits.
 */
class Bits {
public:
    Bits();
    explicit Bits(uint32_t num_bits);

    /*! \brief Set bit at position p. */
    void set(uint32_t p);

    /*! \brief Clear bit at position p. */
    void clear(uint32_t p);

    /*! \brief Toggle bit at position p. */
    void flip(uint32_t p);

    /*! \brief Read bit at position p. */
    bool get(uint32_t p) const;

    /*! \brief Population count (number of set bits). */
    uint32_t count() const;

    /*! \brief Number of bits. */
    uint32_t size() const;

    Bits& operator&=(Bits const& o);
    Bits operator&(Bits const& o) const;

    Bits& operator^=(Bits const& o);
    Bits operator^(Bits const& o) const;
    bool operator==(Bits const& o) const;

    /*! \brief Binary string representation "010101...". */
    std::string to_string() const;

    /*! \brief Raw word access for hashing. */
    std::vector<uint64_t> const& raw() const;

    /*! \brief Mutable raw word access (for fast bulk load). */
    std::vector<uint64_t>& raw_mut();

private:
    uint32_t n_{0};
    std::vector<uint64_t> data_;
};

}  // namespace exact_t

namespace std {
template <>
struct hash<exact_t::Bits> {
    size_t operator()(exact_t::Bits const& b) const {
        size_t h = 14695981039346656037ULL;
        for (uint64_t w : b.raw()) {
            h ^= w;
            h *= 1099511628211ULL;
        }
        return h;
    }
};
}  // namespace std
