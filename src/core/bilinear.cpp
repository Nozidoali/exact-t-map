#include "core/bilinear.hpp"

#include <cassert>
#include <map>
#include <queue>
#include <sstream>
#include <unordered_set>

namespace exact_t {

BilinearFunction::BilinearFunction(uint32_t ni, uint32_t no)
    : num_inputs_(ni), num_outputs_(no) {}

uint32_t BilinearFunction::num_inputs() const { return num_inputs_; }
uint32_t BilinearFunction::num_outputs() const { return num_outputs_; }
std::vector<Triplet> const& BilinearFunction::triplets() const { return triplets_; }

void BilinearFunction::add_triplet(uint32_t i, uint32_t j, uint32_t k) {
    assert(i < num_inputs_);
    assert(j < num_inputs_);
    assert(k < num_outputs_);
    triplets_.push_back({static_cast<uint8_t>(i),
                        static_cast<uint8_t>(j),
                        static_cast<uint8_t>(k)});
}

uint32_t BilinearFunction::total_qubits() const {
    return num_inputs_ + num_outputs_;
}

std::string BilinearFunction::to_string() const {
    std::ostringstream os;
    os << "BilinearFunction(" << num_inputs_ << " -> " << num_outputs_
       << ", " << triplets_.size() << " terms)";
    for (auto const& t : triplets_) {
        os << "\n  x" << static_cast<int>(t.i) << " AND x"
           << static_cast<int>(t.j) << " -> y"
           << static_cast<int>(t.k);
    }
    return os.str();
}

std::vector<uint32_t> BilinearFunction::get_output_parities(uint32_t output_idx) const {
    std::vector<uint32_t> parities;
    for (auto const& t : triplets_) {
        if (t.k == output_idx && t.i < num_inputs_ && t.j < num_inputs_)
            parities.push_back((1u << t.i) | (1u << t.j));
    }
    return parities;
}

namespace {

Bits term_to_vector(uint32_t term, uint32_t vector_dim) {
    Bits vec(vector_dim);
    for (uint32_t t = 1; t < vector_dim; ++t) {
        if ((t & term) == t)
            vec.set(t);
    }
    return vec;
}

}  // anonymous namespace

std::vector<Bits> BilinearFunction::get_kernels(uint32_t dim) const {
    std::vector<Bits> result;
    if (num_outputs_ == 0) return result;

    uint32_t vector_dim = 1u << dim;

    for (uint32_t output_idx = 0; output_idx < num_outputs_; ++output_idx) {
        std::vector<uint32_t> initial_terms = get_output_parities(output_idx);

        if (initial_terms.empty()) {
            result.push_back(Bits(vector_dim));
            continue;
        }

        std::unordered_set<uint32_t> terms_set(initial_terms.begin(), initial_terms.end());
        std::queue<uint32_t> q;
        for (uint32_t term : initial_terms)
            q.push(term);

        while (!q.empty()) {
            uint32_t current = q.front();
            q.pop();
            for (uint32_t other : initial_terms) {
                uint32_t combined = current | other;
                if (combined != 0 && terms_set.find(combined) == terms_set.end()) {
                    terms_set.insert(combined);
                    q.push(combined);
                }
            }
        }

        Bits kernel(vector_dim);
        for (uint32_t term : terms_set)
            kernel ^= term_to_vector(term, vector_dim);
        result.push_back(kernel);
    }

    return result;
}

std::vector<Z8Monomial> BilinearFunction::get_dontcare_monomials(uint32_t dim) const {
    std::vector<Z8Monomial> result;
    if (num_outputs_ == 0) return result;

    for (uint32_t output_idx = 0; output_idx < num_outputs_; ++output_idx) {
        std::vector<uint32_t> monomials = get_output_parities(output_idx);
        Z8Monomial dc_mono(dim);

        if (monomials.empty()) {
            result.push_back(dc_mono);
            continue;
        }

        uint32_t num_mono = static_cast<uint32_t>(monomials.size());
        for (uint32_t mask = 1; mask < (1u << num_mono); ++mask) {
            uint32_t combined_parity = 0;
            uint32_t subset_size = 0;

            for (uint32_t i = 0; i < num_mono; ++i) {
                if (mask & (1u << i)) {
                    combined_parity |= monomials[i];
                    subset_size++;
                }
            }

            int sign = (subset_size % 2 == 1) ? 1 : -1;
            int coeff = sign * (1 << (subset_size - 1));
            dc_mono.add_coeff(combined_parity, coeff);
        }

        result.push_back(dc_mono);
    }

    return result;
}

Tensor BilinearFunction::to_tensor(uint32_t dim) const {
    Tensor t(dim);
    for (auto const& trip : triplets_) {
        assert(trip.i < num_inputs_);
        assert(trip.j < num_inputs_);
        assert(trip.k < num_outputs_);
        t.add_cubic(trip.i, trip.j, num_inputs_ + trip.k);
    }
    return t;
}

Tensor BilinearFunction::to_tensor() const {
    return to_tensor(total_qubits());
}

Phase BilinearFunction::to_phase() const {
    uint32_t nq = total_qubits();
    Phase p(nq);
    for (auto const& t : triplets_)
        p.ccz(t.i, t.j, num_inputs_ + t.k);
    return p;
}

Circuit BilinearFunction::to_circuit() const {
    Circuit circ;
    circ.set_num_qubits(total_qubits());
    for (auto const& t : triplets_) {
        circ.toffoli(t.i, t.j, num_inputs_ + t.k);
    }
    return circ;
}

Circuit BilinearFunction::to_phase_circuit() const {
    uint32_t nq = total_qubits();

    std::map<uint32_t, int> coeffs;
    for (auto const& t : triplets_) {
        uint32_t a = t.i, b = t.j, c = num_inputs_ + t.k;
        uint32_t ma = 1u << a, mb = 1u << b, mc = 1u << c;
        coeffs[ma] = (coeffs[ma] + 1) % 8;
        coeffs[mb] = (coeffs[mb] + 1) % 8;
        coeffs[mc] = (coeffs[mc] + 1) % 8;
        coeffs[ma | mb] = (coeffs[ma | mb] + 7) % 8;
        coeffs[ma | mc] = (coeffs[ma | mc] + 7) % 8;
        coeffs[mb | mc] = (coeffs[mb | mc] + 7) % 8;
        coeffs[ma | mb | mc] = (coeffs[ma | mb | mc] + 1) % 8;
    }

    Circuit circ;
    circ.set_num_qubits(nq);

    for (uint32_t k = 0; k < num_outputs_; ++k)
        circ.h(num_inputs_ + k);

    for (auto const& [mask, coeff] : coeffs) {
        int c = ((coeff % 8) + 8) % 8;
        if (c == 0) continue;

        std::vector<uint32_t> bits;
        for (uint32_t i = 0; i < nq; ++i) {
            if (mask & (1u << i)) bits.push_back(i);
        }
        if (bits.empty()) continue;

        uint32_t target = bits[0];
        for (size_t i = 1; i < bits.size(); ++i)
            circ.cx(bits[i], target);

        circ.append_phase_gates(target, c);

        for (size_t i = bits.size() - 1; i >= 1; --i)
            circ.cx(bits[i], target);
    }

    for (uint32_t k = 0; k < num_outputs_; ++k)
        circ.h(num_inputs_ + k);

    return circ;
}

Tensor to_tensor_with_ancilla_dc(BilinearFunction const& bf, uint32_t dim) {
    assert(dim >= bf.total_qubits());
    Tensor t = bf.to_tensor(dim);
    for (uint32_t q = bf.total_qubits(); q < dim; ++q)
        t.set_output_dont_care(q);
    return t;
}

}  // namespace exact_t
