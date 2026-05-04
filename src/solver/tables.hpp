/*! \file tables.hpp
 *  \brief Precomputed tables for Clifford+T synthesis.
 */

#pragma once

#include <cstdint>
#include <vector>

#include "core/bits.hpp"
#include "core/tensor.hpp"
#include "solver/solver.hpp"

namespace exact_t {

/*! \brief Precomputed lookup tables for parity-to-monomial and Gaussian transforms. */
class SolverTables {
public:
    explicit SolverTables(uint32_t n);

    std::vector<Tensor> const& parity_to_monomial() const;
    std::vector<Z8Phase> const& monomial_to_phase() const;
    Transform const& gaussian_transform() const;
    std::vector<Bits> const& basis_kernels() const;
    uint32_t num_parities() const;

    /*! \brief Save tables to a binary file. Returns true on success. */
    bool save(std::string const& path) const;

    /*! \brief Load tables from a binary file. Returns nullptr if file is
     *  missing, has bad magic/version, or n mismatch. */
    static std::shared_ptr<SolverTables> load(std::string const& path,
                                                uint32_t expected_n);

private:
    SolverTables() = default;  // for deserialization

    uint32_t n_;
    std::vector<Tensor> parity_to_monomial_;
    std::vector<Z8Phase> monomial_to_phase_;
    Transform transform_;
    std::vector<Bits> basis_;
};

}  // namespace exact_t
