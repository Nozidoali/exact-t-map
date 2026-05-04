#include "solver/tables.hpp"

#include <cassert>
#include <cstring>
#include <fstream>

namespace exact_t {

namespace {

constexpr uint32_t TABLES_MAGIC = 0x43544253;  // "STBC" little-endian
constexpr uint32_t TABLES_VERSION = 1;

template <typename T>
void write_pod(std::ostream& os, T const& v) {
    os.write(reinterpret_cast<char const*>(&v), sizeof(T));
}
template <typename T>
bool read_pod(std::istream& is, T& v) {
    is.read(reinterpret_cast<char*>(&v), sizeof(T));
    return is.good();
}

void write_bits(std::ostream& os, Bits const& b) {
    uint32_t n = b.size();
    write_pod(os, n);
    auto const& raw = b.raw();
    uint32_t words = static_cast<uint32_t>(raw.size());
    write_pod(os, words);
    if (words) os.write(reinterpret_cast<char const*>(raw.data()), words * sizeof(uint64_t));
}
bool read_bits(std::istream& is, Bits& b) {
    uint32_t n, words;
    if (!read_pod(is, n) || !read_pod(is, words)) return false;
    b = Bits(n);
    auto& raw = b.raw_mut();
    raw.assign(words, 0);
    if (words) is.read(reinterpret_cast<char*>(raw.data()), words * sizeof(uint64_t));
    return is.good();
}

void write_tensor(std::ostream& os, Tensor const& t) {
    uint32_t n = t.num_vars();
    write_pod(os, n);
    write_bits(os, t.data());
    write_bits(os, t.care());
}
bool read_tensor(std::istream& is, Tensor& t) {
    uint32_t n;
    if (!read_pod(is, n)) return false;
    t = Tensor(n);
    return read_bits(is, t.data()) && read_bits(is, t.care());
}

void write_z8(std::ostream& os, Z8Phase const& z) {
    write_pod(os, z.n);
    uint32_t count = static_cast<uint32_t>(z.terms.size());
    write_pod(os, count);
    for (auto const& [k, v] : z.terms) {
        write_pod(os, k);
        int32_t val = static_cast<int32_t>(v);
        write_pod(os, val);
    }
}
bool read_z8(std::istream& is, Z8Phase& z) {
    uint32_t n, count;
    if (!read_pod(is, n) || !read_pod(is, count)) return false;
    z = Z8Phase(n);
    for (uint32_t i = 0; i < count; ++i) {
        uint32_t k;
        int32_t v;
        if (!read_pod(is, k) || !read_pod(is, v)) return false;
        z.terms[k] = static_cast<int>(v);
    }
    return true;
}

void write_transform(std::ostream& os, Transform const& t) {
    uint32_t mat_count = static_cast<uint32_t>(t.matrix.size());
    write_pod(os, mat_count);
    for (auto const& b : t.matrix) write_bits(os, b);
    uint32_t row_count = static_cast<uint32_t>(t.row_to_cand.size());
    write_pod(os, row_count);
    for (int v : t.row_to_cand) {
        int32_t vi = static_cast<int32_t>(v);
        write_pod(os, vi);
    }
    write_pod(os, t.rank);
}
bool read_transform(std::istream& is, Transform& t) {
    uint32_t mat_count;
    if (!read_pod(is, mat_count)) return false;
    t.matrix.resize(mat_count);
    for (auto& b : t.matrix) if (!read_bits(is, b)) return false;
    uint32_t row_count;
    if (!read_pod(is, row_count)) return false;
    t.row_to_cand.resize(row_count);
    for (auto& v : t.row_to_cand) {
        int32_t vi;
        if (!read_pod(is, vi)) return false;
        v = static_cast<int>(vi);
    }
    return read_pod(is, t.rank);
}

}  // anonymous namespace

namespace {

std::vector<Tensor> gen_parity_to_monomial(uint32_t n) {
    std::vector<Tensor> result;
    result.push_back(Tensor(n));
    for (uint32_t term = 1; term < (1u << n); ++term)
        result.push_back(Tensor::from_parity_mask(n, term));
    return result;
}

std::vector<Z8Phase> gen_monomial_to_phase(uint32_t n) {
    size_t nm = 1u << n;
    std::vector<Z8Phase> result(nm, Z8Phase(n));

    for (size_t m = 1; m < nm; ++m) {
        uint32_t popcount = static_cast<uint32_t>(__builtin_popcount(static_cast<uint32_t>(m)));
        size_t num_subsets = 1u << popcount;

        for (size_t s = 1; s < num_subsets; ++s) {
            uint32_t parity = 0;
            uint32_t bit_idx = 0;
            for (uint32_t i = 0; i < n; ++i) {
                if (m & (1u << i)) {
                    if (s & (1u << bit_idx++))
                        parity |= (1u << i);
                }
            }
            int sign = (__builtin_popcount(static_cast<uint32_t>(s)) % 2) ? 1 : -1;
            result[m].set_coeff(parity, result[m].get_coeff(parity) + sign);
        }
    }
    return result;
}

Transform precompute_gaussian_transform(std::vector<Tensor> const& ptm) {
    Transform transform;
    uint32_t m = static_cast<uint32_t>(ptm.size());
    if (m <= 1) return transform;

    uint32_t nbits = ptm[1].data().size();

    std::vector<Bits> mat;
    std::vector<Bits> bit_rows;
    mat.reserve(nbits);
    bit_rows.reserve(nbits);

    for (uint32_t b = 0; b < nbits; ++b) {
        Bits row(m);
        Bits bit_mask(nbits);
        bit_mask.set(b);
        for (uint32_t parity = 1; parity < m; ++parity) {
            if (ptm[parity].data().get(b))
                row.set(parity);
        }
        mat.push_back(row);
        bit_rows.push_back(bit_mask);
    }

    std::vector<int> pivot_col(nbits, -1);
    uint32_t rank = 0;

    for (uint32_t col = 0; col < m && rank < nbits; ++col) {
        uint32_t pivot_row = rank;
        while (pivot_row < nbits && !mat[pivot_row].get(col))
            ++pivot_row;
        if (pivot_row == nbits) continue;

        if (pivot_row != rank) {
            std::swap(mat[rank], mat[pivot_row]);
            std::swap(bit_rows[rank], bit_rows[pivot_row]);
        }
        pivot_col[rank] = static_cast<int>(col);

        for (uint32_t r = 0; r < nbits; ++r) {
            if (r != rank && mat[r].get(col)) {
                mat[r] ^= mat[rank];
                bit_rows[r] ^= bit_rows[rank];
            }
        }
        ++rank;
    }

    transform.rank = rank;
    for (uint32_t r = 0; r < rank; ++r) {
        transform.matrix.push_back(bit_rows[r]);
        transform.row_to_cand.push_back(pivot_col[r]);
    }
    return transform;
}

std::vector<Bits> compute_basis_kernels(std::vector<Tensor> const& ptm) {
    uint32_t m = static_cast<uint32_t>(ptm.size());
    if (m <= 1) return {};

    uint32_t nbits = ptm[1].data().size();
    std::vector<Bits> mat;
    mat.reserve(nbits);

    for (uint32_t b = 0; b < nbits; ++b) {
        Bits row(m);
        for (uint32_t parity = 1; parity < m; ++parity) {
            if (ptm[parity].data().get(b))
                row.set(parity);
        }
        mat.push_back(row);
    }

    std::vector<int> pivot_col(nbits, -1);
    Bits is_pivot(m);
    uint32_t rank = 0;

    for (uint32_t col = 0; col < m && rank < nbits; ++col) {
        uint32_t pivot_row = rank;
        while (pivot_row < nbits && !mat[pivot_row].get(col))
            ++pivot_row;
        if (pivot_row == nbits) continue;

        if (pivot_row != rank)
            std::swap(mat[rank], mat[pivot_row]);

        pivot_col[rank] = static_cast<int>(col);
        is_pivot.set(col);

        for (uint32_t r = 0; r < nbits; ++r) {
            if (r != rank && mat[r].get(col))
                mat[r] ^= mat[rank];
        }
        ++rank;
    }

    std::vector<Bits> kernels;
    for (uint32_t c = 0; c < m; ++c) {
        if (is_pivot.get(c)) continue;
        Bits k_vec(m);
        k_vec.set(c);
        for (uint32_t r = 0; r < rank; ++r) {
            if (pivot_col[r] >= 0 && mat[r].get(c))
                k_vec.set(static_cast<uint32_t>(pivot_col[r]));
        }
        kernels.push_back(k_vec);
    }
    return kernels;
}

}  // anonymous namespace

SolverTables::SolverTables(uint32_t n)
    : n_(n),
      parity_to_monomial_(gen_parity_to_monomial(n)),
      monomial_to_phase_(gen_monomial_to_phase(n)),
      transform_(precompute_gaussian_transform(parity_to_monomial_)),
      basis_(compute_basis_kernels(parity_to_monomial_)) {
    assert(n <= 16);
}

std::vector<Tensor> const& SolverTables::parity_to_monomial() const {
    return parity_to_monomial_;
}

std::vector<Z8Phase> const& SolverTables::monomial_to_phase() const {
    return monomial_to_phase_;
}

Transform const& SolverTables::gaussian_transform() const {
    return transform_;
}

std::vector<Bits> const& SolverTables::basis_kernels() const {
    return basis_;
}

uint32_t SolverTables::num_parities() const {
    return static_cast<uint32_t>(parity_to_monomial_.size());
}

bool SolverTables::save(std::string const& path) const {
    std::ofstream os(path, std::ios::binary);
    if (!os) return false;
    write_pod(os, TABLES_MAGIC);
    write_pod(os, TABLES_VERSION);
    write_pod(os, n_);
    uint32_t pm_count = static_cast<uint32_t>(parity_to_monomial_.size());
    write_pod(os, pm_count);
    for (auto const& t : parity_to_monomial_) write_tensor(os, t);
    uint32_t mp_count = static_cast<uint32_t>(monomial_to_phase_.size());
    write_pod(os, mp_count);
    for (auto const& z : monomial_to_phase_) write_z8(os, z);
    write_transform(os, transform_);
    uint32_t b_count = static_cast<uint32_t>(basis_.size());
    write_pod(os, b_count);
    for (auto const& b : basis_) write_bits(os, b);
    return os.good();
}

std::shared_ptr<SolverTables> SolverTables::load(std::string const& path,
                                                   uint32_t expected_n) {
    std::ifstream is(path, std::ios::binary);
    if (!is) return nullptr;
    uint32_t magic, version, n;
    if (!read_pod(is, magic) || magic != TABLES_MAGIC) return nullptr;
    if (!read_pod(is, version) || version != TABLES_VERSION) return nullptr;
    if (!read_pod(is, n) || n != expected_n) return nullptr;

    auto tables = std::shared_ptr<SolverTables>(new SolverTables());
    tables->n_ = n;
    uint32_t pm_count;
    if (!read_pod(is, pm_count)) return nullptr;
    tables->parity_to_monomial_.resize(pm_count);
    for (auto& t : tables->parity_to_monomial_) if (!read_tensor(is, t)) return nullptr;
    uint32_t mp_count;
    if (!read_pod(is, mp_count)) return nullptr;
    tables->monomial_to_phase_.resize(mp_count);
    for (auto& z : tables->monomial_to_phase_) if (!read_z8(is, z)) return nullptr;
    if (!read_transform(is, tables->transform_)) return nullptr;
    uint32_t b_count;
    if (!read_pod(is, b_count)) return nullptr;
    tables->basis_.resize(b_count);
    for (auto& b : tables->basis_) if (!read_bits(is, b)) return nullptr;
    return tables;
}

}  // namespace exact_t
