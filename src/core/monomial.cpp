#include "core/monomial.hpp"

#include <sstream>

namespace exact_t {

Z8Monomial::Z8Monomial() = default;

Z8Monomial::Z8Monomial(uint32_t n) : n_(n), coeffs_(1u << n, 0) {}

uint32_t Z8Monomial::num_vars() const { return n_; }

std::vector<int> const& Z8Monomial::coeffs() const { return coeffs_; }

void Z8Monomial::add_coeff(uint32_t idx, int val) {
    coeffs_[idx] = mod8(coeffs_[idx] + val);
}

bool Z8Monomial::is_even() const {
    for (size_t i = 0; i < coeffs_.size(); ++i) {
        if (to_signed(coeffs_[i]) % 2 != 0)
            return false;
    }
    return true;
}

Z8Monomial Z8Monomial::operator+(Z8Monomial const& o) const {
    Z8Monomial r(n_);
    for (uint32_t i = 0; i < static_cast<uint32_t>(coeffs_.size()); ++i)
        r.coeffs_[i] = mod8(coeffs_[i] + o.coeffs_[i]);
    return r;
}

Z8Monomial Z8Monomial::operator-(Z8Monomial const& o) const {
    Z8Monomial r(n_);
    for (uint32_t i = 0; i < static_cast<uint32_t>(coeffs_.size()); ++i)
        r.coeffs_[i] = mod8(coeffs_[i] - o.coeffs_[i]);
    return r;
}

Z8Monomial Z8Monomial::operator*(int scalar) const {
    Z8Monomial r(n_);
    for (uint32_t i = 0; i < static_cast<uint32_t>(coeffs_.size()); ++i)
        r.coeffs_[i] = mod8(coeffs_[i] * scalar);
    return r;
}

std::string Z8Monomial::to_string() const {
    std::ostringstream os;
    bool first = true;
    for (size_t idx = 0; idx < coeffs_.size(); ++idx) {
        if (coeffs_[idx] == 0) continue;
        if (!first) os << "\n";
        first = false;
        bool first_var = true;
        if (idx == 0) {
            os << "  1";
        } else {
            os << "  ";
            for (uint32_t q = 0; q < n_; ++q) {
                if (idx & (1u << q)) {
                    if (!first_var) os << "*";
                    first_var = false;
                    os << "x" << q;
                }
            }
        }
        os << ": " << coeffs_[idx];
    }
    return first ? "  (empty)" : os.str();
}

}  // namespace exact_t
