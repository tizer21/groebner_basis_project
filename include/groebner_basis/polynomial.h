#pragma once

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <map>
#include <ostream>
#include <utility>
#include <vector>

#include "monomial.h"
#include "monomial_orders.h"
#include "term.h"

namespace groebner_basis {

template <typename Coefficient, typename MonomialOrder = LexicographicOrder>
class Polynomial {
public:
    Polynomial() = default;

    Polynomial(std::initializer_list<Term<Coefficient>> terms) : terms_(BuildTermContainer(terms)) {
        assert(IsInvariantSatisfied());
    }

    explicit Polynomial(const std::vector<Term<Coefficient>>& terms)
        : terms_(BuildTermContainer(terms)) {
        assert(IsInvariantSatisfied());
    }

    std::vector<Term<Coefficient>> GetTerms() const {
        std::vector<Term<Coefficient>> result;
        result.reserve(terms_.size());
        for (const auto& [monomial, coefficient] : terms_) {
            result.emplace_back(coefficient, monomial);
        }
        return result;
    }

    const Coefficient& GetLeadingCoefficient() const {
        assert(!IsZero() && "The zero polynomial has no leading coefficient");
        return terms_.begin()->second;
    }

    const Monomial& GetLeadingMonomial() const {
        assert(!IsZero() && "The zero polynomial has no leading monomial");
        return terms_.begin()->first;
    }

    Term<Coefficient> GetLeadingTerm() const {
        assert(!IsZero() && "The zero polynomial has no leading term");
        return Term<Coefficient>{GetLeadingCoefficient(), GetLeadingMonomial()};
    }

    bool IsZero() const {
        return terms_.empty();
    }

    bool IsOne() const {
        return terms_.size() == 1 && GetLeadingMonomial().IsOne() &&
               GetLeadingCoefficient() == Coefficient{1};
    }

    bool IsMonic() const {
        return !IsZero() && GetLeadingCoefficient() == Coefficient{1};
    }

    Polynomial& MakeMonic();

    Polynomial operator-() const {
        Polynomial result = *this;
        for (auto& [monomial, coefficient] : result.terms_) {
            coefficient = -coefficient;
        }

        assert(result.IsInvariantSatisfied());
        return result;
    }

    Polynomial& operator+=(const Coefficient& addend);
    Polynomial& operator+=(const Term<Coefficient>& addend);
    Polynomial& operator+=(const Polynomial& addend);

    Polynomial& operator-=(const Coefficient& subtrahend);
    Polynomial& operator-=(const Term<Coefficient>& subtrahend);
    Polynomial& operator-=(const Polynomial& subtrahend);

    Polynomial& operator*=(const Coefficient& factor);
    Polynomial& operator*=(const Term<Coefficient>& factor);
    Polynomial& operator*=(const Polynomial& factor);

    friend Polynomial operator+(Polynomial lhs, const Coefficient& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend Polynomial operator+(const Coefficient& lhs, Polynomial rhs) {
        rhs += lhs;
        return rhs;
    }

    friend Polynomial operator+(Polynomial lhs, const Term<Coefficient>& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend Polynomial operator+(const Term<Coefficient>& lhs, Polynomial rhs) {
        rhs += lhs;
        return rhs;
    }

    friend Polynomial operator+(Polynomial lhs, const Polynomial& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend Polynomial operator-(Polynomial lhs, const Coefficient& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend Polynomial operator-(const Coefficient& lhs, Polynomial rhs) {
        rhs *= -Coefficient{1};
        rhs += lhs;
        return rhs;
    }

    friend Polynomial operator-(Polynomial lhs, const Term<Coefficient>& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend Polynomial operator-(const Term<Coefficient>& lhs, Polynomial rhs) {
        rhs *= -Coefficient{1};
        rhs += lhs;
        return rhs;
    }

    friend Polynomial operator-(Polynomial lhs, const Polynomial& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend Polynomial operator*(Polynomial lhs, const Coefficient& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Polynomial operator*(const Coefficient& lhs, Polynomial rhs) {
        rhs *= lhs;
        return rhs;
    }

    friend Polynomial operator*(Polynomial lhs, const Term<Coefficient>& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Polynomial operator*(const Term<Coefficient>& lhs, Polynomial rhs) {
        rhs *= lhs;
        return rhs;
    }

    friend Polynomial operator*(Polynomial lhs, const Polynomial& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend bool operator==(const Polynomial& lhs, const Polynomial& rhs) = default;

    friend std::ostream& operator<<(std::ostream& os, const Polynomial& polynomial) {
        return polynomial.WriteTo(os);
    }

private:
    using TermMap = std::map<Monomial, Coefficient, MonomialOrder>;

    static void AddTermTo(
        TermMap& terms, const Monomial& monomial, const Coefficient& coefficient
    ) {
        if (coefficient == Coefficient{0}) {
            return;
        }

        auto [position, is_inserted] = terms.try_emplace(monomial, coefficient);
        if (!is_inserted) {
            position->second += coefficient;
        }
        if (position->second == Coefficient{0}) {
            terms.erase(position);
        }
    }

    template <typename TermRange>
    static TermMap BuildTermContainer(const TermRange& terms) {
        TermMap result;
        for (const auto& term : terms) {
            AddTermTo(result, term.GetMonomial(), term.GetCoefficient());
        }

        return result;
    }

    bool IsInvariantSatisfied() const {
        return std::all_of(terms_.begin(), terms_.end(), [](const auto& entry) {
            return entry.second != Coefficient{0};
        });
    }

    void AddTerm(const Monomial& monomial, const Coefficient& coefficient) {
        AddTermTo(terms_, monomial, coefficient);
    }

    std::ostream& WriteTo(std::ostream& os) const;

    TermMap terms_;
};

// Implementation details only below here.

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder>& Polynomial<Coefficient, MonomialOrder>::MakeMonic() {
    if (IsZero() || IsMonic()) {
        return *this;
    }

    const auto leading_coefficient = GetLeadingCoefficient();
    for (auto& [monomial, coefficient] : terms_) {
        coefficient /= leading_coefficient;
    }

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder>& Polynomial<Coefficient, MonomialOrder>::operator+=(
    const Coefficient& addend
) {
    AddTerm({}, addend);

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder>& Polynomial<Coefficient, MonomialOrder>::operator+=(
    const Term<Coefficient>& addend
) {
    AddTerm(addend.GetMonomial(), addend.GetCoefficient());

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder>& Polynomial<Coefficient, MonomialOrder>::operator+=(
    const Polynomial& addend
) {
    if (this == &addend) {
        const Polynomial copy = addend;
        return *this += copy;
    }

    for (const auto& [monomial, coefficient] : addend.terms_) {
        AddTerm(monomial, coefficient);
    }

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder>& Polynomial<Coefficient, MonomialOrder>::operator-=(
    const Coefficient& subtrahend
) {
    AddTerm({}, -subtrahend);

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder>& Polynomial<Coefficient, MonomialOrder>::operator-=(
    const Term<Coefficient>& subtrahend
) {
    AddTerm(subtrahend.GetMonomial(), -subtrahend.GetCoefficient());

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder>& Polynomial<Coefficient, MonomialOrder>::operator-=(
    const Polynomial& subtrahend
) {
    if (this == &subtrahend) {
        terms_.clear();

        assert(IsInvariantSatisfied());
        return *this;
    }

    for (const auto& [monomial, coefficient] : subtrahend.terms_) {
        AddTerm(monomial, -coefficient);
    }

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder>& Polynomial<Coefficient, MonomialOrder>::operator*=(
    const Coefficient& factor
) {
    const Coefficient factor_copy = factor;

    if (factor_copy == Coefficient{0}) {
        terms_.clear();

        assert(IsInvariantSatisfied());
        return *this;
    }

    if (factor_copy == Coefficient{1}) {
        return *this;
    }

    for (auto& [monomial, coefficient] : terms_) {
        coefficient *= factor_copy;
    }

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder>& Polynomial<Coefficient, MonomialOrder>::operator*=(
    const Term<Coefficient>& factor
) {
    if (IsZero()) {
        return *this;
    }

    if (factor.IsZero()) {
        terms_.clear();
        assert(IsInvariantSatisfied());
        return *this;
    }

    if (factor.GetMonomial().IsOne()) {
        *this *= factor.GetCoefficient();
        return *this;
    }

    TermMap result;
    for (const auto& [monomial, coefficient] : terms_) {
        result.emplace_hint(
            result.end(), monomial * factor.GetMonomial(), coefficient * factor.GetCoefficient()
        );
    }
    terms_ = std::move(result);

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder>& Polynomial<Coefficient, MonomialOrder>::operator*=(
    const Polynomial& factor
) {
    if (IsZero()) {
        return *this;
    }

    if (factor.IsZero()) {
        terms_.clear();

        assert(IsInvariantSatisfied());
        return *this;
    }

    if (factor.IsOne()) {
        return *this;
    }

    if (IsOne()) {
        *this = factor;
        return *this;
    }

    if (factor.terms_.size() == 1) {
        const Term<Coefficient> factor_term = factor.GetLeadingTerm();
        return *this *= factor_term;
    }

    Polynomial result;
    for (const auto& [monomial, coefficient] : terms_) {
        for (const auto& [other_monomial, other_coefficient] : factor.terms_) {
            result.AddTerm(monomial * other_monomial, coefficient * other_coefficient);
        }
    }
    *this = std::move(result);

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient, typename MonomialOrder>
std::ostream& Polynomial<Coefficient, MonomialOrder>::WriteTo(std::ostream& os) const {
    if (IsZero()) {
        return os << '0';
    }

    bool is_first = true;
    for (const auto& [monomial, coefficient] : terms_) {
        if (is_first) {
            os << Term<Coefficient>{coefficient, monomial};

            is_first = false;
            continue;
        }

        if (coefficient < Coefficient{0}) {
            os << " - " << Term<Coefficient>{-coefficient, monomial};
        } else {
            os << " + " << Term<Coefficient>{coefficient, monomial};
        }
    }

    return os;
}

}  // namespace groebner_basis
