#pragma once

#include <stddef.h>

#include <algorithm>
#include <cassert>
#include <set>
#include <utility>
#include <vector>

#include "polynomial.h"
#include "term.h"

namespace groebner_basis::algorithms {

template <typename Coefficient, typename MonomialOrder>
using PolynomialBasis = std::vector<Polynomial<Coefficient, MonomialOrder>>;

template <typename Coefficient, typename MonomialOrder>
struct PolynomialDivisionResult {
    std::vector<Polynomial<Coefficient, MonomialOrder>> quotients;
    Polynomial<Coefficient, MonomialOrder> remainder;
};

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder> ComputeSPolynomial(
    const Polynomial<Coefficient, MonomialOrder>& lhs,
    const Polynomial<Coefficient, MonomialOrder>& rhs
);

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder> ReducePolynomial(
    Polynomial<Coefficient, MonomialOrder> polynomial,
    const std::vector<Polynomial<Coefficient, MonomialOrder>>& divisors
);

namespace detail {

template <bool kComputeQuotients, typename Coefficient, typename MonomialOrder>
PolynomialDivisionResult<Coefficient, MonomialOrder> DividePolynomialImpl(
    Polynomial<Coefficient, MonomialOrder> dividend,
    const std::vector<Polynomial<Coefficient, MonomialOrder>>& divisors
) {
    PolynomialDivisionResult<Coefficient, MonomialOrder> result;
    if constexpr (kComputeQuotients) {
        result.quotients.resize(divisors.size());
    }

    while (!dividend.IsZero()) {
        const auto dividend_lt = dividend.GetLeadingTerm();

        bool was_reduced = false;
        for (size_t idx = 0; idx != divisors.size(); ++idx) {
            const auto& divisor = divisors[idx];
            if (divisor.IsZero()) {
                continue;
            }

            const auto divisor_lt = divisor.GetLeadingTerm();
            if (!divisor_lt.Divides(dividend_lt)) {
                continue;
            }

            const auto quotient_term = dividend_lt / divisor_lt;
            if constexpr (kComputeQuotients) {
                result.quotients[idx] += quotient_term;
            }

            dividend -= quotient_term * divisor;

            was_reduced = true;
            break;
        }

        if (!was_reduced) {
            result.remainder += dividend_lt;
            dividend -= dividend_lt;
        }
    }

    return result;
}

using BuchbergerPair = std::pair<size_t, size_t>;
using BuchbergerPairSet = std::set<BuchbergerPair>;

inline BuchbergerPair MakeBuchbergerPair(size_t lhs, size_t rhs) {
    assert(lhs != rhs);

    if (lhs < rhs) {
        return {lhs, rhs};
    }
    return {rhs, lhs};
}

template <typename Coefficient, typename MonomialOrder>
BuchbergerPairSet BuildBuchbergerPairs(const PolynomialBasis<Coefficient, MonomialOrder>& basis) {
    BuchbergerPairSet result;

    for (size_t lhs_idx = 0; lhs_idx != basis.size(); ++lhs_idx) {
        for (size_t rhs_idx = lhs_idx + 1; rhs_idx != basis.size(); ++rhs_idx) {
            result.emplace_hint(result.end(), lhs_idx, rhs_idx);
        }
    }

    return result;
}

template <typename Coefficient, typename MonomialOrder>
bool MeetsBuchbergerProductCriterion(
    const Polynomial<Coefficient, MonomialOrder>& lhs,
    const Polynomial<Coefficient, MonomialOrder>& rhs
) {
    const auto& lhs_leading_monomial = lhs.GetLeadingMonomial();
    const auto& rhs_leading_monomial = rhs.GetLeadingMonomial();

    const auto least_common_multiple =
        ComputeLeastCommonMultiple(lhs_leading_monomial, rhs_leading_monomial);

    return least_common_multiple == lhs_leading_monomial * rhs_leading_monomial;
}

template <typename Coefficient, typename MonomialOrder>
bool MeetsBuchbergerChainCriterion(
    size_t lhs_idx, size_t rhs_idx, const PolynomialBasis<Coefficient, MonomialOrder>& basis,
    const BuchbergerPairSet& pending_pairs
) {
    const auto least_common_multiple = ComputeLeastCommonMultiple(
        basis[lhs_idx].GetLeadingMonomial(), basis[rhs_idx].GetLeadingMonomial()
    );

    for (size_t middle_idx = 0; middle_idx != basis.size(); ++middle_idx) {
        if (middle_idx == lhs_idx || middle_idx == rhs_idx) {
            continue;
        }

        const auto lhs_middle_pair = MakeBuchbergerPair(lhs_idx, middle_idx);
        const auto rhs_middle_pair = MakeBuchbergerPair(rhs_idx, middle_idx);

        if (pending_pairs.contains(lhs_middle_pair) || pending_pairs.contains(rhs_middle_pair)) {
            continue;
        }

        if (basis[middle_idx].GetLeadingMonomial().Divides(least_common_multiple)) {
            return true;
        }
    }

    return false;
}

enum class Improvement : unsigned char {
    kNone,
    kProductCriterion,
    kProductAndChainCriteria,
};

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> RunBuchbergerAlgorithm(
    std::vector<Polynomial<Coefficient, MonomialOrder>> generators,
    Improvement improvement = Improvement::kProductAndChainCriteria
) {
    auto basis = std::move(generators);
    basis.erase(
        std::remove_if(
            basis.begin(), basis.end(), [](const auto& polynomial) { return polynomial.IsZero(); }
        ),
        basis.end()
    );

    auto pending_pairs = BuildBuchbergerPairs(basis);

    while (!pending_pairs.empty()) {
        const auto pair_position = pending_pairs.begin();
        const auto [lhs_idx, rhs_idx] = *pair_position;
        pending_pairs.erase(pair_position);

        if (improvement != Improvement::kNone &&
            MeetsBuchbergerProductCriterion(basis[lhs_idx], basis[rhs_idx])) {
            continue;
        }

        if (improvement == Improvement::kProductAndChainCriteria &&
            MeetsBuchbergerChainCriterion(lhs_idx, rhs_idx, basis, pending_pairs)) {
            continue;
        }

        auto remainder =
            ReducePolynomial(ComputeSPolynomial(basis[lhs_idx], basis[rhs_idx]), basis);

        if (remainder.IsZero()) {
            continue;
        }

        const size_t new_idx = basis.size();
        basis.push_back(std::move(remainder));
        for (size_t existing_idx = 0; existing_idx != new_idx; ++existing_idx) {
            pending_pairs.emplace(existing_idx, new_idx);
        }
    }

    return basis;
}

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> MinimizeGroebnerBasis(
    PolynomialBasis<Coefficient, MonomialOrder> basis
) {
    std::vector<bool> should_keep(basis.size(), true);
    for (size_t candidate_idx = 0; candidate_idx != basis.size(); ++candidate_idx) {
        const auto& candidate_lm = basis[candidate_idx].GetLeadingMonomial();

        for (size_t current_idx = 0; current_idx != basis.size(); ++current_idx) {
            if (candidate_idx == current_idx) {
                continue;
            }

            const auto& current_lm = basis[current_idx].GetLeadingMonomial();

            const bool is_candidate_redundant =
                current_lm.Divides(candidate_lm) &&
                (current_lm != candidate_lm || current_idx < candidate_idx);
            if (is_candidate_redundant) {
                should_keep[candidate_idx] = false;
                break;
            }
        }
    }

    size_t write_idx = 0;
    for (size_t read_idx = 0; read_idx != basis.size(); ++read_idx) {
        if (!should_keep[read_idx]) {
            continue;
        }

        if (write_idx != read_idx) {
            basis[write_idx] = std::move(basis[read_idx]);
        }

        ++write_idx;
    }
    basis.resize(write_idx);

    for (auto& polynomial : basis) {
        polynomial.MakeMonic();

        if (polynomial.IsOne()) {
            return {std::move(polynomial)};
        }
    }

    return basis;
}

template <typename Coefficient, typename MonomialOrder>
void SortBasisByLeadingMonomial(PolynomialBasis<Coefficient, MonomialOrder>& basis) {
    std::stable_sort(basis.begin(), basis.end(), [](const auto& lhs, const auto& rhs) {
        return MonomialOrder{}(lhs.GetLeadingMonomial(), rhs.GetLeadingMonomial());
    });
}

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> InterreduceGroebnerBasis(
    PolynomialBasis<Coefficient, MonomialOrder> minimal_basis
) {
    for (size_t idx = 0; idx != minimal_basis.size(); ++idx) {
        auto dividend = std::move(minimal_basis[idx]);
        minimal_basis[idx] = Polynomial<Coefficient, MonomialOrder>{};

        auto remainder = ReducePolynomial(std::move(dividend), minimal_basis);

        assert(!remainder.IsZero());
        assert(remainder.IsMonic());

        minimal_basis[idx] = std::move(remainder);
    }

    SortBasisByLeadingMonomial(minimal_basis);
    return minimal_basis;
}

}  // namespace detail

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder> ComputeSPolynomial(
    const Polynomial<Coefficient, MonomialOrder>& lhs,
    const Polynomial<Coefficient, MonomialOrder>& rhs
) {
    assert(!lhs.IsZero() && !rhs.IsZero() && "S-polynomial requires two nonzero polynomials");

    const auto least_common_multiple =
        ComputeLeastCommonMultiple(lhs.GetLeadingMonomial(), rhs.GetLeadingMonomial());
    const Term<Coefficient> common_term{least_common_multiple};

    const auto lhs_factor = common_term / lhs.GetLeadingTerm();
    const auto rhs_factor = common_term / rhs.GetLeadingTerm();

    return lhs * lhs_factor - rhs * rhs_factor;
}

template <typename Coefficient, typename MonomialOrder>
PolynomialDivisionResult<Coefficient, MonomialOrder> DividePolynomial(
    Polynomial<Coefficient, MonomialOrder> dividend,
    const std::vector<Polynomial<Coefficient, MonomialOrder>>& divisors
) {
    return detail::DividePolynomialImpl<true>(std::move(dividend), divisors);
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder> ReducePolynomial(
    Polynomial<Coefficient, MonomialOrder> polynomial,
    const std::vector<Polynomial<Coefficient, MonomialOrder>>& divisors
) {
    auto result = detail::DividePolynomialImpl<false>(std::move(polynomial), divisors);
    return std::move(result.remainder);
}

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> ComputeGroebnerBasis(
    std::vector<Polynomial<Coefficient, MonomialOrder>> generators
) {
    return detail::RunBuchbergerAlgorithm(std::move(generators));
}

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> ComputeMinimalGroebnerBasis(
    std::vector<Polynomial<Coefficient, MonomialOrder>> generators
) {
    return detail::MinimizeGroebnerBasis(ComputeGroebnerBasis(std::move(generators)));
}

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> ComputeReducedGroebnerBasis(
    std::vector<Polynomial<Coefficient, MonomialOrder>> generators
) {
    return detail::InterreduceGroebnerBasis(ComputeMinimalGroebnerBasis(std::move(generators)));
}

}  // namespace groebner_basis::algorithms
