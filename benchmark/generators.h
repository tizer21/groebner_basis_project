#pragma once

#include <stddef.h>

#include <cassert>
#include <utility>
#include <vector>

#include "groebner_basis/monomial.h"
#include "groebner_basis/polynomial.h"
#include "groebner_basis/term.h"

namespace groebner_basis::benchmarking {

template <typename Coefficient, typename MonomialOrder>
std::vector<Polynomial<Coefficient, MonomialOrder>> GenerateCyclicSystem(size_t variable_count) {
    assert(variable_count > 0);

    using PolynomialType = Polynomial<Coefficient, MonomialOrder>;
    using TermType = Term<Coefficient>;

    std::vector<PolynomialType> cyclic_system;
    cyclic_system.reserve(variable_count);

    for (size_t product_length = 1; product_length < variable_count; ++product_length) {
        std::vector<TermType> terms;
        terms.reserve(variable_count);

        for (size_t start_idx = 0; start_idx < variable_count; ++start_idx) {
            std::vector<Monomial::Degree> degrees(variable_count, Monomial::Degree{0});
            for (size_t offset = 0; offset < product_length; ++offset) {
                const size_t variable_idx = (start_idx + offset) % variable_count;
                degrees[variable_idx] = Monomial::Degree{1};
            }

            terms.emplace_back(Coefficient{1}, Monomial{std::move(degrees)});
        }

        cyclic_system.emplace_back(terms);
    }

    std::vector<TermType> final_terms;
    final_terms.reserve(2);

    std::vector<Monomial::Degree> final_degrees(variable_count, Monomial::Degree{1});
    final_terms.emplace_back(Coefficient{1}, Monomial{std::move(final_degrees)});
    final_terms.emplace_back(Coefficient{-1});

    cyclic_system.emplace_back(final_terms);

    return cyclic_system;
}

template <typename Coefficient, typename MonomialOrder>
std::vector<Polynomial<Coefficient, MonomialOrder>> GenerateKatsuraSystem(size_t order) {
    assert(order > 0);

    using PolynomialType = Polynomial<Coefficient, MonomialOrder>;
    using TermType = Term<Coefficient>;

    const size_t variable_count = order + 1;

    std::vector<PolynomialType> katsura_system;
    katsura_system.reserve(variable_count);

    std::vector<TermType> normalization_terms;
    normalization_terms.reserve(variable_count + 1);

    for (size_t variable_idx = 0; variable_idx < variable_count; ++variable_idx) {
        std::vector<Monomial::Degree> degrees(variable_count, Monomial::Degree{0});
        degrees[variable_idx] = Monomial::Degree{1};

        const Coefficient coefficient = variable_idx == 0 ? Coefficient{1} : Coefficient{2};

        normalization_terms.emplace_back(coefficient, Monomial{std::move(degrees)});
    }
    normalization_terms.emplace_back(Coefficient{-1});
    katsura_system.emplace_back(normalization_terms);

    for (size_t equation_idx = 0; equation_idx < order; ++equation_idx) {
        std::vector<TermType> terms;
        terms.reserve(2 * variable_count);

        const auto addQuadraticTerm = [&](size_t lhs_idx, size_t rhs_idx) {
            std::vector<Monomial::Degree> degrees(variable_count, Monomial::Degree{0});

            ++degrees[lhs_idx];
            ++degrees[rhs_idx];

            terms.emplace_back(Coefficient{1}, Monomial{std::move(degrees)});
        };

        for (size_t lhs_idx = 0; lhs_idx < variable_count; ++lhs_idx) {
            const size_t rhs_idx =
                lhs_idx <= equation_idx ? equation_idx - lhs_idx : lhs_idx - equation_idx;
            addQuadraticTerm(lhs_idx, rhs_idx);
        }

        for (size_t reflected_idx = 1; reflected_idx <= order - equation_idx; ++reflected_idx) {
            addQuadraticTerm(reflected_idx, equation_idx + reflected_idx);
        }

        std::vector<Monomial::Degree> subtrahend_degrees(variable_count, Monomial::Degree{0});
        subtrahend_degrees[equation_idx] = Monomial::Degree{1};
        terms.emplace_back(Coefficient{-1}, Monomial{std::move(subtrahend_degrees)});

        katsura_system.emplace_back(terms);
    }

    return katsura_system;
}

}  // namespace groebner_basis::benchmarking
