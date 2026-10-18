#pragma once

#include <cassert>
#include <compare>
#include <cstdint>
#include <numeric>
#include <ostream>

namespace groebner_basis::fields {

class Rational {
public:
    using Integer = int64_t;

    Rational() = default;

    explicit Rational(Integer numerator) : Rational(numerator, Integer{1}) {}

    explicit Rational(Integer numerator, Integer denominator);

    constexpr Integer GetNumerator() const noexcept {
        return numerator_;
    }

    constexpr Integer GetDenominator() const noexcept {
        return denominator_;
    }

    constexpr bool IsZero() const noexcept {
        return numerator_ == 0;
    }

    Rational operator-() const {
        return Rational{-numerator_, denominator_};
    }

    Rational& operator+=(const Rational& addend);
    Rational& operator-=(const Rational& subtrahend);
    Rational& operator*=(const Rational& factor);
    Rational& operator/=(const Rational& divisor);

    friend Rational operator+(Rational lhs, const Rational& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend Rational operator-(Rational lhs, const Rational& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend Rational operator*(Rational lhs, const Rational& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Rational operator/(Rational lhs, const Rational& rhs) {
        lhs /= rhs;
        return lhs;
    }

    friend constexpr bool operator==(const Rational& lhs, const Rational& rhs) noexcept = default;

    friend constexpr std::strong_ordering operator<=>(
        const Rational& lhs, const Rational& rhs
    ) noexcept {
        const Integer lhs_scaled = lhs.numerator_ * rhs.denominator_;
        const Integer rhs_scaled = rhs.numerator_ * lhs.denominator_;

        return lhs_scaled <=> rhs_scaled;
    }

private:
    struct CanonicalRepresentation {
        Integer numerator;
        Integer denominator;
    };

    explicit Rational(CanonicalRepresentation representation);

    bool IsInvariantSatisfied() const {
        return denominator_ > 0 && std::gcd(numerator_, denominator_) == Integer{1};
    }

    static CanonicalRepresentation Canonicalize(Integer numerator, Integer denominator);

    Integer numerator_{0};
    Integer denominator_{1};
};

inline std::ostream& operator<<(std::ostream& out, const Rational& rational);

// Implementation details only below here.

inline Rational::Rational(Integer numerator, Integer denominator)
    : Rational(Canonicalize(numerator, denominator)) {
    assert(IsInvariantSatisfied());
}

inline Rational::Rational(CanonicalRepresentation representation)
    : numerator_(representation.numerator), denominator_(representation.denominator) {
    assert(IsInvariantSatisfied());
}

inline Rational::CanonicalRepresentation Rational::Canonicalize(
    Integer numerator, Integer denominator
) {
    assert(denominator != 0 && "The denominator must not be zero");

    const Integer divisor = std::gcd(numerator, denominator);
    numerator /= divisor;
    denominator /= divisor;

    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }

    return {numerator, denominator};
}

inline Rational& Rational::operator+=(const Rational& addend) {
    const auto new_representation = Canonicalize(
        numerator_ * addend.denominator_ + addend.numerator_ * denominator_,
        denominator_ * addend.denominator_
    );

    numerator_ = new_representation.numerator;
    denominator_ = new_representation.denominator;

    assert(IsInvariantSatisfied());
    return *this;
}

inline Rational& Rational::operator-=(const Rational& subtrahend) {
    const auto new_representation = Canonicalize(
        numerator_ * subtrahend.denominator_ - subtrahend.numerator_ * denominator_,
        denominator_ * subtrahend.denominator_
    );

    numerator_ = new_representation.numerator;
    denominator_ = new_representation.denominator;

    assert(IsInvariantSatisfied());
    return *this;
}

inline Rational& Rational::operator*=(const Rational& factor) {
    const auto new_representation =
        Canonicalize(numerator_ * factor.numerator_, denominator_ * factor.denominator_);

    numerator_ = new_representation.numerator;
    denominator_ = new_representation.denominator;

    assert(IsInvariantSatisfied());
    return *this;
}

inline Rational& Rational::operator/=(const Rational& divisor) {
    assert(!divisor.IsZero() && "The divisor must not be zero");

    const auto new_representation =
        Canonicalize(numerator_ * divisor.denominator_, denominator_ * divisor.numerator_);

    numerator_ = new_representation.numerator;
    denominator_ = new_representation.denominator;

    assert(IsInvariantSatisfied());
    return *this;
}

inline std::ostream& operator<<(std::ostream& out, const Rational& rational) {
    out << rational.GetNumerator();

    if (rational.GetDenominator() != 1) {
        out << '/' << rational.GetDenominator();
    }

    return out;
}

}  // namespace groebner_basis::fields
