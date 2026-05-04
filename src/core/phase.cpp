#include "core/phase.hpp"

#include <algorithm>
#include <cassert>
#include <sstream>
#include <vector>

namespace exact_t {

Phase::Phase() = default;

Phase::Phase(uint32_t n) : n_(n) {}

uint32_t Phase::num_vars() const { return n_; }

std::set<uint32_t> const& Phase::terms() const { return terms_; }

void Phase::toggle_term(uint32_t t) {
    t &= (1u << n_) - 1;
    std::set<uint32_t>::iterator it = terms_.find(t);
    if (it != terms_.end())
        terms_.erase(it);
    else
        terms_.insert(t);
}

void Phase::ccz(uint32_t i, uint32_t j, uint32_t k) {
    assert(i < n_ && j < n_ && k < n_);
    uint32_t mi = 1u << i, mj = 1u << j, mk = 1u << k;
    toggle_term(mi | mj | mk);
    toggle_term(mi | mj);
    toggle_term(mi | mk);
    toggle_term(mj | mk);
    toggle_term(mi);
    toggle_term(mj);
    toggle_term(mk);
}

uint32_t Phase::size() const { return static_cast<uint32_t>(terms_.size()); }

bool Phase::empty() const { return terms_.empty(); }

bool Phase::contains(uint32_t t) const { return terms_.count(t) > 0; }

void Phase::apply_inverse_cnot(uint32_t ctrl, uint32_t target) {
    assert(ctrl < n_ && target < n_ && ctrl != target);
    std::set<uint32_t> new_terms;
    for (uint32_t t : terms_) {
        uint32_t nt = t;
        if (t & (1u << target))
            nt ^= (1u << ctrl);
        new_terms.insert(nt);
    }
    terms_ = new_terms;
}

Phase& Phase::operator+=(Phase const& o) {
    for (uint32_t t : o.terms_)
        toggle_term(t);
    return *this;
}

uint64_t Phase::hash() const {
    if (terms_.empty()) return 0;

    uint32_t num_terms = static_cast<uint32_t>(terms_.size());
    std::vector<uint32_t> term_vec(terms_.begin(), terms_.end());

    std::vector<std::vector<bool>> columns(n_, std::vector<bool>(num_terms, false));
    for (uint32_t t = 0; t < num_terms; ++t) {
        for (uint32_t v = 0; v < n_; ++v) {
            if (term_vec[t] & (1u << v))
                columns[v][t] = true;
        }
    }

    std::sort(columns.begin(), columns.end());

    std::vector<uint32_t> canonical(num_terms, 0);
    for (uint32_t v = 0; v < n_; ++v) {
        for (uint32_t t = 0; t < num_terms; ++t) {
            if (columns[v][t])
                canonical[t] |= (1u << v);
        }
    }
    std::sort(canonical.begin(), canonical.end());

    uint64_t h = 14695981039346656037ULL;
    for (uint32_t c : canonical) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    return h;
}

std::string Phase::to_string() const {
    if (terms_.empty()) return "0";
    std::ostringstream os;
    bool first = true;
    for (uint32_t mask : terms_) {
        if (!first) os << " + ";
        first = false;
        if (mask == 0) {
            os << "1";
            continue;
        }
        bool first_var = true;
        for (uint32_t i = 0; i < n_; ++i) {
            if (mask & (1u << i)) {
                if (!first_var) os << "*";
                first_var = false;
                os << "z" << i;
            }
        }
    }
    return os.str();
}

}  // namespace exact_t
