#include "solver/solver.hpp"
#include "solver/gauss.hpp"
#include "solver/tables.hpp"
#include "core/bilinear.hpp"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>

namespace exact_t {

Z8Phase::Z8Phase() = default;
Z8Phase::Z8Phase(uint32_t n) : n(n) {}

void Z8Phase::set_coeff(uint32_t parity, int coeff) {
    coeff = Z8Monomial::mod8(coeff);
    if (coeff == 0)
        terms.erase(parity);
    else
        terms[parity] = coeff;
}

int Z8Phase::get_coeff(uint32_t parity) const {
    std::unordered_map<uint32_t, int>::const_iterator it = terms.find(parity);
    return it != terms.end() ? it->second : 0;
}

void Z8Phase::add_coeff(uint32_t parity, int val) {
    set_coeff(parity, get_coeff(parity) + val);
}

bool Z8Phase::empty() const { return terms.empty(); }

uint32_t Z8Phase::size() const { return static_cast<uint32_t>(terms.size()); }

std::vector<std::pair<uint32_t, int>> Z8Phase::to_vector() const {
    std::vector<std::pair<uint32_t, int>> result;
    for (auto const& [parity, coeff] : terms)
        result.push_back({parity, coeff});
    std::sort(result.begin(), result.end());
    return result;
}

Z8Phase Z8Phase::from_phase(Phase const& phase) {
    Z8Phase result(phase.num_vars());
    for (uint32_t parity : phase.terms())
        result.add_coeff(parity, 1);
    return result;
}

Z8Phase Z8Phase::operator+(Z8Phase const& o) const {
    Z8Phase result(n);
    for (auto const& [parity, coeff] : terms)
        result.set_coeff(parity, coeff);
    for (auto const& [parity, coeff] : o.terms)
        result.add_coeff(parity, coeff);
    return result;
}

std::pair<Z8Phase, Z8Phase> Z8Phase::filter_by_qubits(uint32_t qubit_mask) const {
    Z8Phase involved(n), not_involved(n);
    for (auto const& [parity, coeff] : terms) {
        if (parity == 0 || (parity & qubit_mask) != 0)
            involved.set_coeff(parity, coeff);
        else
            not_involved.set_coeff(parity, coeff);
    }
    return {involved, not_involved};
}

namespace {

uint32_t count_odd_z8_coefficients(Z8Phase const& phase) {
    uint32_t count = 0;
    for (auto const& [parity, coeff] : phase.terms) {
        if (coeff % 2 != 0)
            count++;
    }
    return count;
}

uint32_t clean_ancilla_compensation_t_count(std::vector<int> const& dc_coeffs) {
    uint32_t count = 0;
    for (int c : dc_coeffs) {
        int comp = ((8 - c) % 8 + 8) % 8;
        if (comp % 2 != 0)
            count++;
    }
    return count;
}

std::vector<std::string> split_csv(std::string const& line) {
    std::vector<std::string> result;
    std::istringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ','))
        result.push_back(token);
    return result;
}

Tensor tensor_from_bitstring(std::string const& s, uint32_t n) {
    Tensor t(n);
    for (size_t i = 0; i < s.length() && i < t.data().size(); ++i)
        if (s[i] == '1') t.data().flip(static_cast<uint32_t>(i));
    return t;
}

Phase phase_from_solution_bits(std::string const& bits, uint32_t n) {
    Phase p(n);
    for (size_t parity = 1; parity < bits.size(); ++parity)
        if (bits[parity] == '1')
            p.toggle_term(static_cast<uint32_t>(parity));
    return p;
}

void dump_n7_instance(Tensor const& target, uint32_t heuristic_popcount) {
    static char const* const path = std::getenv("EXACT_T_DUMP_N7_INSTANCES");
    if (path == nullptr) return;
    if (target.num_vars() != 7) return;

    static std::ofstream ofs(path, std::ios::app);
    if (!ofs.is_open()) return;

    ofs << target.data().to_string() << ',';
    if (!target.has_dont_cares()) {
        ofs << std::string(target.data().size(), '1');
    } else {
        ofs << target.care().to_string();
    }
    ofs << ',' << heuristic_popcount << '\n';
    ofs.flush();
}

}  // anonymous namespace

CliffordTSolver::CliffordTSolver(uint32_t n, SolverConfig const& config)
    : max_n_(n), config_(config) {
    assert(n <= 16);
}

std::shared_ptr<SolverTables> CliffordTSolver::get_tables(uint32_t n) const {
    auto it = tables_cache_.find(n);
    if (it != tables_cache_.end()) return it->second;

    std::shared_ptr<SolverTables> tables;
    char const* cache_dir = std::getenv("EXACT_T_TABLES_CACHE_DIR");
    std::string cache_path;
    if (cache_dir != nullptr) {
        cache_path = std::string(cache_dir) + "/tables_n" + std::to_string(n) + ".bin";
        tables = SolverTables::load(cache_path, n);
    }
    if (!tables) {
        tables = std::make_shared<SolverTables>(n);
        if (!cache_path.empty()) tables->save(cache_path);
    }
    tables_cache_[n] = tables;
    return tables;
}

void CliffordTSolver::load_cache(std::string const& path) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) return;

    std::string line;
    if (!std::getline(ifs, line)) return;
    auto meta = split_csv(line);
    if (meta.size() < 2 || meta[0] != "n") return;
    uint32_t file_n = static_cast<uint32_t>(std::strtoul(meta[1].c_str(), nullptr, 10));

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        auto cols = split_csv(line);
        if (cols.size() < 4) continue;

        Tensor tensor = tensor_from_bitstring(cols[0], file_n);
        Phase phase = phase_from_solution_bits(cols[3], file_n);
        cache_[tensor] = phase;
    }
}

Z8Monomial CliffordTSolver::phase_to_z8(Phase const& phase, uint32_t n) const {
    auto tables = get_tables(n);
    Z8Monomial result(n);
    for (uint32_t parity : phase.terms()) {
        Z8Monomial parity_mono = tables->parity_to_monomial()[parity].to_z8_monomial();
        for (size_t monomial = 0; monomial < (1u << n); ++monomial)
            result.add_coeff(static_cast<uint32_t>(monomial), parity_mono[static_cast<uint32_t>(monomial)]);
    }
    return result;
}

Z8Phase CliffordTSolver::solve_z8(Z8Monomial const& diff, uint32_t n) const {
    auto tables = get_tables(n);
    size_t nm = 1u << n;
    Z8Phase z8(n);

    for (size_t p = 0; p < nm; ++p) {
        int x_p = 0;
        for (size_t m = 0; m < nm; ++m) {
            int coeff = tables->monomial_to_phase()[m].get_coeff(static_cast<uint32_t>(p));
            if (coeff != 0) {
                int divisor = Z8Monomial::degree_divisor(static_cast<uint32_t>(m));
                int b_m = Z8Monomial::to_signed(diff[static_cast<uint32_t>(m)]);
                x_p = (x_p + coeff * b_m / divisor) % 8;
            }
        }
        if (x_p != 0)
            z8.set_coeff(static_cast<uint32_t>(p), x_p);
    }
    return z8;
}

std::vector<std::pair<uint32_t, uint32_t>> CliffordTSolver::greedy_cnot_preprocess(
    Tensor& target, uint32_t n) const {
    std::vector<std::pair<uint32_t, uint32_t>> applied_cnots;

    bool improved_cnot = true;
    while (improved_cnot) {
        improved_cnot = false;
        int best_ctrl = -1, best_tgt = -1;
        uint32_t best_terms = target.num_terms();

        for (uint32_t ctrl = 0; ctrl < n; ++ctrl) {
            for (uint32_t tgt = 0; tgt < n; ++tgt) {
                if (ctrl == tgt) continue;
                Tensor transformed = target;
                transformed.apply_cnot(ctrl, tgt);
                uint32_t nt = transformed.num_terms();
                if (nt < best_terms) {
                    best_terms = nt;
                    best_ctrl = static_cast<int>(ctrl);
                    best_tgt = static_cast<int>(tgt);
                }
            }
        }

        if (best_ctrl >= 0) {
            target.apply_cnot(static_cast<uint32_t>(best_ctrl),
                              static_cast<uint32_t>(best_tgt));
            applied_cnots.push_back({static_cast<uint32_t>(best_ctrl),
                                     static_cast<uint32_t>(best_tgt)});
            improved_cnot = true;
        }
    }

    return applied_cnots;
}

void CliffordTSolver::emit_phase_rotations(Circuit& circ, Z8Phase const& z8_phase,
                                             uint32_t n) const {
    for (auto const& [mask, coeff] : z8_phase.terms) {
        int c = Z8Monomial::mod8(coeff);
        if (c == 0) continue;

        std::vector<uint32_t> bits;
        for (uint32_t i = 0; i < n; ++i) {
            if (mask & (1u << i)) bits.push_back(i);
        }
        if (bits.empty()) continue;

        uint32_t tgt = bits[0];
        for (size_t i = 1; i < bits.size(); ++i)
            circ.cx(bits[i], tgt);

        circ.append_phase_gates(tgt, c);

        for (size_t i = bits.size() - 1; i >= 1; --i)
            circ.cx(bits[i], tgt);
    }
}

Phase CliffordTSolver::solve(Tensor const& target) const {
    std::unordered_map<Tensor, Phase>::const_iterator cache_it = cache_.find(target);
    if (cache_it != cache_.end()) {
        cache_hits_++;
        dump_n7_instance(target, cache_it->second.size());
        return cache_it->second;
    }

    uint32_t n = target.num_vars();
    auto tables = get_tables(n);

    Tensor solve_target = target;
    auto applied_cnots = greedy_cnot_preprocess(solve_target, n);

    uint32_t m = tables->num_parities();
    Transform const& transform = tables->gaussian_transform();
    std::vector<Bits> const& basis = tables->basis_kernels();

    auto solve_single = [&](Tensor const& t) -> Phase {
        Tensor t_care = t;
        if (t.has_dont_cares()) {
            uint32_t nbits = t.data().size();
            for (uint32_t b = 0; b < nbits; ++b) {
                if (!t.is_care(b) && t.data().get(b))
                    t_care.data().flip(b);
            }
        }

        std::vector<int> sol = GaussianSolver::solve(t_care, transform);
        std::vector<int> sol_backup = sol;

        std::vector<Bits> augmented_basis = basis;
        auto dc_freedom = GaussianSolver::dc_freedom(t, transform, m);
        for (auto const& fv : dc_freedom)
            augmented_basis.push_back(fv);

        GaussianSolver::refine_basis(augmented_basis, sol, m, config_.max_flip_depth);

        Phase refined = GaussianSolver::to_phase(sol, n);
        if (!GaussianSolver::verify(refined, t))
            return GaussianSolver::to_phase(sol_backup, n);
        return refined;
    };

    Phase phase_orig = solve_single(target);
    Phase phase = phase_orig;

    if (!applied_cnots.empty()) {
        Phase phase_cnot = solve_single(solve_target);
        for (auto it = applied_cnots.rbegin(); it != applied_cnots.rend(); ++it)
            phase_cnot.apply_inverse_cnot(it->first, it->second);

        if (phase_cnot.size() < phase_orig.size())
            phase = phase_cnot;
    }

    assert(GaussianSolver::verify(phase, target));
    cache_misses_++;
    cache_[target] = phase;
    dump_n7_instance(target, phase.size());
    return phase;
}

Circuit CliffordTSolver::synthesize(BilinearFunction const& bf) const {
    uint32_t n = bf.total_qubits();
    assert(n <= max_n_);

    Tensor target = bf.to_tensor(n);
    Phase phase = solve(target);

    Z8Monomial target_z8 = target.to_z8_monomial();
    Z8Monomial synth_z8 = phase_to_z8(phase, n);
    Z8Monomial diff = target_z8 - synth_z8;
    assert(diff.is_even());

    Z8Phase corrections = solve_z8(diff, n);
    Z8Phase z8_phase = Z8Phase::from_phase(phase) + corrections;

    Circuit circ;
    circ.set_num_qubits(n);

    for (uint32_t k = 0; k < bf.num_outputs(); ++k)
        circ.h(bf.num_inputs() + k);

    emit_phase_rotations(circ, z8_phase, n);

    for (uint32_t k = 0; k < bf.num_outputs(); ++k)
        circ.h(bf.num_inputs() + k);

    return circ;
}

Phase CliffordTSolver::solve_with_kernels(Tensor const& target,
                                           std::vector<Bits> const& extra_kernels) const {
    std::unordered_map<Tensor, Phase>::const_iterator cache_it = cache_.find(target);
    if (cache_it != cache_.end()) {
        cache_hits_++;
        return cache_it->second;
    }

    uint32_t n = target.num_vars();
    auto tables = get_tables(n);

    Tensor solve_target = target;
    auto applied_cnots = greedy_cnot_preprocess(solve_target, n);

    uint32_t m = tables->num_parities();
    Transform const& transform = tables->gaussian_transform();
    std::vector<Bits> const& basis = tables->basis_kernels();

    std::vector<Bits> augmented_basis = basis;
    for (auto const& kernel : extra_kernels) {
        Bits candidate_kernel(m);
        for (uint32_t i = 0; i < kernel.size() && i < m; ++i) {
            if (kernel.get(i))
                candidate_kernel.set(i);
        }
        augmented_basis.push_back(candidate_kernel);
    }

    auto solve_single = [&](Tensor const& t) -> Phase {
        std::vector<int> sol = GaussianSolver::solve(t, transform);
        std::vector<int> sol_backup = sol;
        GaussianSolver::refine_basis(augmented_basis, sol, m, config_.max_flip_depth);

        Phase refined = GaussianSolver::to_phase(sol, n);
        if (!GaussianSolver::verify(refined, t))
            return GaussianSolver::to_phase(sol_backup, n);
        return refined;
    };

    Phase phase_orig = solve_single(target);
    Phase phase = phase_orig;

    if (!applied_cnots.empty()) {
        Phase phase_cnot = solve_single(solve_target);
        for (auto it = applied_cnots.rbegin(); it != applied_cnots.rend(); ++it)
            phase_cnot.apply_inverse_cnot(it->first, it->second);

        if (phase_cnot.size() < phase_orig.size())
            phase = phase_cnot;
    }

    assert(GaussianSolver::verify(phase, target));
    cache_misses_++;
    cache_[target] = phase;
    return phase;
}

Circuit CliffordTSolver::synthesize(BilinearFunction const& bf,
                                    bool use_ancilla_dc) const {
    if (!use_ancilla_dc)
        return synthesize(bf);

    uint32_t n = max_n_;
    assert(bf.total_qubits() <= n);

    Tensor target = bf.to_tensor(n);

    std::vector<Bits> kernels = bf.get_kernels(n);
    Phase phase = solve_with_kernels(target, kernels);

    Z8Monomial target_z8 = target.to_z8_monomial();
    Z8Monomial synth_z8 = phase_to_z8(phase, n);
    Z8Monomial diff_z8 = target_z8 - synth_z8;

    std::vector<Z8Monomial> dc_monomials = bf.get_dontcare_monomials(n);

    std::vector<uint32_t> output_qubits;
    for (uint32_t k = 0; k < bf.num_outputs(); ++k)
        output_qubits.push_back(bf.num_inputs() + k);

    uint32_t output_qubit_mask = 0;
    for (uint32_t q : output_qubits)
        output_qubit_mask |= (1u << q);

    Z8Monomial best_diff = diff_z8;
    std::vector<int> best_coeffs(dc_monomials.size(), 0);
    uint32_t best_t = UINT32_MAX;

    uint32_t num_dc = static_cast<uint32_t>(dc_monomials.size());
    if (num_dc <= 3) {
        uint32_t num_combos = 1;
        for (uint32_t i = 0; i < num_dc; ++i)
            num_combos *= 8;

        for (uint32_t combo = 0; combo < num_combos; ++combo) {
            std::vector<int> dc_coeffs(num_dc, 0);
            uint32_t temp = combo;
            for (uint32_t i = 0; i < num_dc; ++i) {
                dc_coeffs[i] = temp % 8;
                temp /= 8;
            }

            Z8Monomial trial_diff = diff_z8;
            for (uint32_t i = 0; i < num_dc; ++i) {
                if (dc_coeffs[i] != 0)
                    trial_diff = trial_diff - (dc_monomials[i] * dc_coeffs[i]);
            }

            if (!trial_diff.is_even())
                continue;

            Z8Phase corrections = solve_z8(trial_diff, n);
            Z8Phase z8_phase = Z8Phase::from_phase(phase) + corrections;
            auto [filtered, non_output] = z8_phase.filter_by_qubits(output_qubit_mask);

            uint32_t t_count = count_odd_z8_coefficients(filtered)
                             + count_odd_z8_coefficients(non_output)
                             + clean_ancilla_compensation_t_count(dc_coeffs);

            if (t_count < best_t) {
                best_t = t_count;
                best_diff = trial_diff;
                best_coeffs = dc_coeffs;
            }
        }
    }

    if (best_t == UINT32_MAX) {
        best_diff = diff_z8;
        best_coeffs.assign(num_dc, 0);
    }

    assert(best_diff.is_even());

    Z8Phase corrections = solve_z8(best_diff, n);
    Z8Phase z8_phase = Z8Phase::from_phase(phase) + corrections;

    auto [filtered_z8, non_output_z8] = z8_phase.filter_by_qubits(output_qubit_mask);

    Circuit circ;
    circ.set_num_qubits(n);

    emit_phase_rotations(circ, non_output_z8, n);

    for (uint32_t k = 0; k < bf.num_outputs(); ++k)
        circ.h(bf.num_inputs() + k);

    emit_phase_rotations(circ, filtered_z8, n);

    for (uint32_t k = 0; k < bf.num_outputs(); ++k)
        circ.h(bf.num_inputs() + k);

    for (uint32_t i = 0; i < num_dc; ++i) {
        int comp = ((8 - best_coeffs[i]) % 8 + 8) % 8;
        if (comp == 0) continue;
        circ.append_phase_gates(output_qubits[i], comp);
    }

    return circ;
}

Circuit CliffordTSolver::synthesize_tensor(Tensor const& target,
                                            uint32_t num_inputs,
                                            uint32_t num_outputs) const {
    uint32_t n = target.num_vars();
    assert(n <= max_n_);

    Phase phase = solve(target);

    Z8Monomial target_z8 = target.to_z8_monomial();
    Z8Monomial synth_z8 = phase_to_z8(phase, n);
    Z8Monomial diff = target_z8 - synth_z8;

    if (target.has_dont_cares()) {
        uint32_t nm = 1u << n;
        for (uint32_t m = 0; m < nm; ++m) {
            if (!target.is_care(m))
                diff[m] = 0;
        }
    }

    assert(diff.is_even());

    Z8Phase corrections = solve_z8(diff, n);
    Z8Phase z8_phase = Z8Phase::from_phase(phase) + corrections;

    Circuit circ;
    circ.set_num_qubits(n);

    for (uint32_t k = 0; k < num_outputs; ++k)
        circ.h(num_inputs + k);

    emit_phase_rotations(circ, z8_phase, n);

    for (uint32_t k = 0; k < num_outputs; ++k)
        circ.h(num_inputs + k);

    return circ;
}

Circuit synthesize_clifford_t(CliffordTSolver const& solver,
                               BilinearFunction const& bf) {
    return solver.synthesize(bf);
}

}  // namespace exact_t
