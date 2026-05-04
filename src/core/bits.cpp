#include "core/bits.hpp"

#include <cassert>

namespace exact_t {

Bits::Bits() = default;

Bits::Bits(uint32_t num_bits) : n_(num_bits), data_((num_bits + 63) / 64, 0) {}

void Bits::set(uint32_t p) {
    assert(p < n_);
    data_[p / 64] |= 1ULL << (p % 64);
}

void Bits::clear(uint32_t p) {
    assert(p < n_);
    data_[p / 64] &= ~(1ULL << (p % 64));
}

void Bits::flip(uint32_t p) {
    assert(p < n_);
    data_[p / 64] ^= 1ULL << (p % 64);
}

bool Bits::get(uint32_t p) const {
    assert(p < n_);
    return (data_[p / 64] >> (p % 64)) & 1;
}

uint32_t Bits::count() const {
    uint32_t c = 0;
    for (uint64_t w : data_)
        c += static_cast<uint32_t>(__builtin_popcountll(w));
    return c;
}

uint32_t Bits::size() const { return n_; }

Bits& Bits::operator&=(Bits const& o) {
    for (size_t i = 0; i < data_.size() && i < o.data_.size(); ++i)
        data_[i] &= o.data_[i];
    return *this;
}

Bits Bits::operator&(Bits const& o) const {
    Bits result = *this;
    result &= o;
    return result;
}

Bits& Bits::operator^=(Bits const& o) {
    assert(n_ == o.n_);
    for (size_t i = 0; i < data_.size(); ++i)
        data_[i] ^= o.data_[i];
    return *this;
}

Bits Bits::operator^(Bits const& o) const {
    Bits r = *this;
    r ^= o;
    return r;
}

bool Bits::operator==(Bits const& o) const {
    return n_ == o.n_ && data_ == o.data_;
}

std::string Bits::to_string() const {
    std::string s(n_, '0');
    for (uint32_t i = 0; i < n_; ++i)
        if (get(i)) s[i] = '1';
    return s;
}

std::vector<uint64_t> const& Bits::raw() const { return data_; }
std::vector<uint64_t>& Bits::raw_mut() { return data_; }

}  // namespace exact_t
