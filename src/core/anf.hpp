/*! \file anf.hpp
 *  \brief Algebraic Normal Form representation up to degree 2.
 *
 *  Represents Boolean functions: f = c + sum(xi) + sum(xi*xj)
 *  \tparam VarType Variable type (must support == and <)
 */

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace exact_t {

/*! \brief Algebraic Normal Form up to degree 2 over GF(2). */
template <typename VarType>
class ANF2 {
    static_assert(
        std::is_same_v<decltype(std::declval<VarType>() == std::declval<VarType>()), bool>,
        "VarType must support operator==");
    static_assert(
        std::is_same_v<decltype(std::declval<VarType>() < std::declval<VarType>()), bool>,
        "VarType must support operator<");

  public:
    bool degree0 = false;
    std::vector<VarType> degree1;
    std::vector<std::pair<VarType, VarType>> degree2;

    ANF2() = default;
    ANF2(bool constant) : degree0(constant) {}
    ANF2(VarType const& var) { add_degree1_term(var); }

    /*! \brief Check if ANF is constant (no variables). */
    bool is_constant() const { return degree1.empty() && degree2.empty(); }

    /*! \brief Get constant value (only valid if is_constant()). */
    bool get_constant_value(bool complement = false) const {
        if (!is_constant()) return false;
        bool val = !degree0;
        return complement ? !val : val;
    }

    /*! \brief Create constant ANF. */
    static ANF2 make_constant(bool value) { return ANF2(value); }

    /*! \brief Clear all terms. */
    void clear() {
        degree0 = false;
        degree1.clear();
        degree2.clear();
    }

    /*! \brief Add linear term. */
    void add_degree1_term(VarType var) { degree1.push_back(var); }

    /*! \brief Add quadratic term. */
    void add_degree2_term(VarType i, VarType j) {
        degree2.emplace_back(i, j);
    }

    /*! \brief Get polynomial degree (0, 1, or 2). */
    uint32_t get_degree() const {
        if (!degree2.empty()) return 2;
        if (!degree1.empty()) return 1;
        return 0;
    }

    /*! \brief Convert to string representation. */
    std::string to_string() const {
        std::string result;
        if (!degree1.empty()) {
            result += "degree1_terms=[";
            for (size_t i = 0; i < degree1.size(); ++i) {
                if (i > 0) result += ", ";
                result += std::to_string(degree1[i]);
            }
            result += "]";
        }
        if (!degree2.empty()) {
            if (!result.empty()) result += ", ";
            result += "degree2_terms=[";
            for (size_t i = 0; i < degree2.size(); ++i) {
                if (i > 0) result += ", ";
                result += "(" + std::to_string(degree2[i].first) + "," +
                          std::to_string(degree2[i].second) + ")";
            }
            result += "]";
        }
        return result;
    }

    /*! \brief XOR two ANFs (addition in GF(2)). */
    ANF2 operator+(ANF2 const& other) const {
        ANF2 result;
        result.degree0 = degree0 ^ other.degree0;

        std::map<VarType, uint32_t> d1_count;
        for (auto const& t : degree1) d1_count[t]++;
        for (auto const& t : other.degree1) d1_count[t]++;
        for (auto const& [t, c] : d1_count)
            if (c % 2 == 1) result.degree1.push_back(t);

        std::map<std::pair<VarType, VarType>, uint32_t> d2_count;
        for (auto const& t : degree2) d2_count[t]++;
        for (auto const& t : other.degree2) d2_count[t]++;
        for (auto const& [t, c] : d2_count)
            if (c % 2 == 1) result.degree2.push_back(t);

        return result;
    }

    /*! \brief AND two ANFs (multiplication in GF(2)). */
    ANF2 operator*(ANF2 const& other) const {
        ANF2 result;
        result.degree0 = degree0 && other.degree0;

        for (auto const& t0 : degree1) {
            for (auto const& t1 : other.degree1) {
                if (t0 == t1) {
                    result.degree0 = !result.degree0;
                } else {
                    auto pair = t0 < t1 ? std::make_pair(t0, t1)
                                        : std::make_pair(t1, t0);
                    result.degree2.push_back(pair);
                }
            }
        }

        if (degree0)
            for (auto const& t : other.degree1) result.degree1.push_back(t);
        if (other.degree0)
            for (auto const& t : degree1) result.degree1.push_back(t);

        std::map<std::pair<VarType, VarType>, uint32_t> d2_count;
        for (auto const& t : result.degree2) d2_count[t]++;
        result.degree2.clear();
        for (auto const& [t, c] : d2_count)
            if (c % 2 == 1) result.degree2.push_back(t);

        std::map<VarType, uint32_t> d1_count;
        for (auto const& t : result.degree1) d1_count[t]++;
        result.degree1.clear();
        for (auto const& [t, c] : d1_count)
            if (c % 2 == 1) result.degree1.push_back(t);

        return result;
    }

    /*! \brief Complement ANF (toggle constant term). */
    ANF2 operator!() const {
        ANF2 result = *this;
        result.degree0 = !result.degree0;
        return result;
    }

    /*! \brief Merge two ANFs with optional complements and operation. */
    static ANF2 merge(ANF2 const& a, ANF2 const& b,
                      bool compl_a, bool compl_b, bool is_and) {
        ANF2 a_compl = compl_a ? !a : a;
        ANF2 b_compl = compl_b ? !b : b;
        return is_and ? a_compl * b_compl : a_compl + b_compl;
    }
};

}  // namespace exact_t
