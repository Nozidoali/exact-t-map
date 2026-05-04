#include "core/tensor.hpp"
#include "core/monomial.hpp"
#include "core/phase.hpp"

#include <cassert>

namespace exact_t {

Tensor::Tensor() = default;

Tensor::Tensor(uint32_t n) : n_(n), data_(1u << n) {
    assert(n <= 20);
}

uint32_t Tensor::num_vars() const { return n_; }

Bits const& Tensor::data() const { return data_; }
Bits& Tensor::data() { return data_; }

Bits const& Tensor::care() const { return care_; }
Bits& Tensor::care() { return care_; }

void Tensor::set_dont_care(uint32_t mask) {
    if (care_.size() == 0) {
        care_ = Bits(1u << n_);
        for (uint32_t i = 0; i < (1u << n_); ++i)
            care_.set(i);
    }
    care_.clear(mask);
}

void Tensor::set_output_dont_care(uint32_t output_bit) {
    assert(output_bit < n_);
    uint32_t num_monomials = 1u << n_;
    for (uint32_t mask = 0; mask < num_monomials; ++mask) {
        if (mask & (1u << output_bit))
            set_dont_care(mask);
    }
}

bool Tensor::is_care(uint32_t mask) const {
    if (care_.size() == 0) return true;
    return care_.get(mask);
}

bool Tensor::has_dont_cares() const {
    if (care_.size() == 0) return false;
    return care_.count() < care_.size();
}

uint32_t Tensor::num_dont_cares() const {
    if (care_.size() == 0) return 0;
    return care_.size() - care_.count();
}

void Tensor::add_linear(uint32_t i) {
    assert(i < n_);
    data_.flip(1u << i);
}

void Tensor::add_quadratic(uint32_t i, uint32_t j) {
    assert(i < n_ && j < n_ && i != j);
    data_.flip((1u << i) | (1u << j));
}

void Tensor::add_cubic(uint32_t i, uint32_t j, uint32_t k) {
    assert(i < n_ && j < n_ && k < n_);
    assert(i != j && i != k && j != k);
    data_.flip((1u << i) | (1u << j) | (1u << k));
}

void Tensor::apply_cnot(uint32_t ctrl, uint32_t target) {
    assert(ctrl < n_ && target < n_ && ctrl != target);
    uint32_t num_monomials = 1u << n_;
    Bits new_data(num_monomials);
    Bits new_care;
    bool has_dc = has_dont_cares();
    if (has_dc)
        new_care = Bits(num_monomials);

    for (uint32_t mask = 0; mask < num_monomials; ++mask) {
        uint32_t new_mask = mask;
        if (mask & (1u << target))
            new_mask ^= (1u << ctrl);

        if (data_.get(mask))
            new_data.flip(new_mask);

        if (has_dc && care_.get(mask))
            new_care.set(new_mask);
    }
    data_ = new_data;
    if (has_dc)
        care_ = new_care;
}

uint32_t Tensor::num_terms() const {
    return data_.count();
}

std::string Tensor::to_string() const {
    return data_.to_string();
}

Tensor& Tensor::operator+=(Tensor const& o) {
    assert(n_ == o.n_);
    data_ ^= o.data_;
    return *this;
}

Tensor Tensor::operator+(Tensor const& o) const {
    Tensor r = *this;
    r += o;
    return r;
}

bool Tensor::operator==(Tensor const& o) const {
    if (n_ != o.n_ || !(data_ == o.data_)) return false;
    bool a_dc = has_dont_cares(), b_dc = o.has_dont_cares();
    if (a_dc != b_dc) return false;
    if (a_dc && !(care_ == o.care_)) return false;
    return true;
}

Tensor Tensor::from_parity_mask(uint32_t n, uint32_t mask) {
    Tensor t(n);
    std::vector<uint32_t> bits;
    for (uint32_t i = 0; i < n; ++i) {
        if (mask & (1u << i)) bits.push_back(i);
    }
    uint32_t deg = static_cast<uint32_t>(bits.size());
    for (uint32_t a = 0; a < deg; ++a) {
        for (uint32_t b = a; b < deg; ++b) {
            for (uint32_t c = b; c < deg; ++c) {
                uint32_t m = (1u << bits[a]) | (1u << bits[b]) | (1u << bits[c]);
                t.data_.set(m);
            }
        }
    }
    return t;
}

Z8Monomial Tensor::to_z8_monomial() const {
    Z8Monomial result(n_);
    uint32_t num_monomials = 1u << n_;
    for (uint32_t mask = 0; mask < num_monomials; ++mask) {
        if (data_.get(mask))
            result[mask] = Z8Monomial::t_coefficient(mask);
    }
    return result;
}

Tensor Tensor::from_phase(Phase const& p) {
    Tensor t(p.num_vars());
    for (uint32_t term : p.terms())
        t += from_parity_mask(p.num_vars(), term);
    return t;
}

}  // namespace exact_t
