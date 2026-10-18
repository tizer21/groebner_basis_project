#pragma once

#include <cassert>
#include <ostream>
#include <utility>

#include "monomial.h"

namespace groebner_basis {

template <typename Coefficient>
class Term {
public:
    Term() = default;
    explicit Term(Coefficient coefficient);
    explicit Term(Monomial monomial);
    explicit Term(Coefficient coefficient, Monomial monomial);

    const Coefficient& GetCoefficient() const noexcept {
        return coefficient_;
    }

    const Monomial& GetMonomial() const noexcept {
        return monomial_;
    }

    bool IsZero() const {
        return coefficient_ == Coefficient{0};
    }

    bool Divides(const Term& dividend) const {
        if (IsZero()) {
            return false;
        }
        if (dividend.IsZero()) {
            return true;
        }
        return monomial_.Divides(dividend.monomial_);
    }

    Term operator-() const {
        return Term{-coefficient_, monomial_};
    }

    Term& operator*=(const Coefficient& factor);
    Term& operator*=(const Monomial& factor);
    Term& operator*=(const Term& factor);
    Term& operator/=(const Coefficient& divisor);
    Term& operator/=(const Monomial& divisor);
    Term& operator/=(const Term& divisor);

    friend Term operator*(Term lhs, const Coefficient& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Term operator*(const Coefficient& lhs, Term rhs) {
        rhs *= lhs;
        return rhs;
    }

    friend Term operator*(Term lhs, const Monomial& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Term operator*(const Monomial& lhs, Term rhs) {
        rhs *= lhs;
        return rhs;
    }

    friend Term operator*(Term lhs, const Term& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Term operator/(Term lhs, const Coefficient& rhs) {
        lhs /= rhs;
        return lhs;
    }

    friend Term operator/(Term lhs, const Monomial& rhs) {
        lhs /= rhs;
        return lhs;
    }

    friend Term operator/(Term lhs, const Term& rhs) {
        lhs /= rhs;
        return lhs;
    }

    friend bool operator==(const Term&, const Term&) = default;

private:
    static Monomial CanonicalizeMonomialForCoefficient(
        Monomial&& monomial, const Coefficient& coefficient
    ) {
        if (coefficient == Coefficient{0}) {
            return Monomial{};
        }
        return monomial;
    }

    bool IsInvariantSatisfied() const {
        return !IsZero() || monomial_.IsOne();
    }

    Coefficient coefficient_{0};
    Monomial monomial_;
};

template <typename Coefficient>
std::ostream& operator<<(std::ostream& os, const Term<Coefficient>& term);

// Implementation details only below here.

template <typename Coefficient>
Term<Coefficient>::Term(Coefficient coefficient) : Term(std::move(coefficient), Monomial{}) {}

template <typename Coefficient>
Term<Coefficient>::Term(Monomial monomial) : Term(Coefficient{1}, std::move(monomial)) {}

template <typename Coefficient>
Term<Coefficient>::Term(Coefficient coefficient, Monomial monomial)
    : coefficient_(std::move(coefficient)),
      monomial_(CanonicalizeMonomialForCoefficient(std::move(monomial), coefficient_)) {
    assert(IsInvariantSatisfied());
}

template <typename Coefficient>
Term<Coefficient>& Term<Coefficient>::operator*=(const Coefficient& factor) {
    if (IsZero()) {
        return *this;
    }

    if (factor == Coefficient{0}) {
        coefficient_ = Coefficient{0};
        monomial_ = Monomial{};
    } else if (factor != Coefficient{1}) {
        coefficient_ *= factor;
    }

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient>
Term<Coefficient>& Term<Coefficient>::operator*=(const Monomial& factor) {
    if (!IsZero() && !factor.IsOne()) {
        monomial_ *= factor;
    }

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient>
Term<Coefficient>& Term<Coefficient>::operator*=(const Term& factor) {
    if (IsZero()) {
        return *this;
    }

    if (factor.IsZero()) {
        coefficient_ = Coefficient{0};
        monomial_ = Monomial{};

        assert(IsInvariantSatisfied());
        return *this;
    }

    if (factor.coefficient_ != Coefficient{1}) {
        coefficient_ *= factor.coefficient_;
    }

    if (IsZero()) {
        monomial_ = Monomial{};
    } else if (!factor.monomial_.IsOne()) {
        monomial_ *= factor.monomial_;
    }

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient>
Term<Coefficient>& Term<Coefficient>::operator/=(const Coefficient& divisor) {
    assert(divisor != Coefficient{0} && "Term division requires a nonzero coefficient");

    coefficient_ /= divisor;

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient>
Term<Coefficient>& Term<Coefficient>::operator/=(const Monomial& divisor) {
    if (!IsZero()) {
        assert(
            divisor.Divides(monomial_) &&
            "Term division requires that divisor monomial divides dividend"
        );

        monomial_ /= divisor;
    }

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient>
Term<Coefficient>& Term<Coefficient>::operator/=(const Term& divisor) {
    assert(divisor.Divides(*this) && "Term division requires that divisor divides dividend");

    if (!IsZero()) {
        coefficient_ /= divisor.coefficient_;
        monomial_ /= divisor.monomial_;
    }

    assert(IsInvariantSatisfied());
    return *this;
}

template <typename Coefficient>
std::ostream& operator<<(std::ostream& os, const Term<Coefficient>& term) {
    const auto& coefficient = term.GetCoefficient();
    const auto& monomial = term.GetMonomial();

    if (term.IsZero()) {
        return os << coefficient;
    }

    if (monomial.IsOne()) {
        return os << coefficient;
    }

    const Coefficient one{1};
    if (coefficient == one) {
        return os << monomial;
    }

    const Coefficient zero{0};
    if (coefficient < zero && coefficient == -one) {
        return os << '-' << monomial;
    }

    return os << coefficient << '*' << monomial;
}

}  // namespace groebner_basis
