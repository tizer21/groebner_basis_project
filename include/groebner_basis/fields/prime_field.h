#pragma once

#include <cassert>
#include <compare>
#include <cstdint>
#include <ostream>

namespace groebner_basis::fields {
namespace detail {

constexpr bool IsPrime(int64_t number) {
    if (number < 2) {
        return false;
    }

    for (int64_t i = 2; i * i <= number; ++i) {
        if (number % i == 0) {
            return false;
        }
    }

    return true;
}

inline int64_t MultiplyModulo(int64_t lhs, int64_t rhs, int64_t modulus) {
    return lhs * rhs % modulus;
}

inline int64_t ComputeModularExponentiation(int64_t base, int64_t exponent, int64_t modulus) {
    base %= modulus;
    if (base < 0) {
        base += modulus;
    }

    int64_t result = 1;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = MultiplyModulo(result, base, modulus);
        }
        base = MultiplyModulo(base, base, modulus);
        exponent /= 2;
    }

    return result;
}

}  // namespace detail

template <int64_t modulus>
class PrimeField {
    static_assert(modulus > 1, "Modulus must be greater than 1");
    static_assert(detail::IsPrime(modulus), "Modulus must be prime");

public:
    using Integer = int64_t;

    PrimeField() = default;

    explicit PrimeField(Integer value);

    constexpr Integer GetValue() const noexcept {
        return value_;
    }

    PrimeField operator-() const;

    PrimeField& operator+=(const PrimeField& addend);
    PrimeField& operator-=(const PrimeField& subtrahend);
    PrimeField& operator*=(const PrimeField& factor);
    PrimeField& operator/=(const PrimeField& divisor);

    friend PrimeField operator+(PrimeField lhs, const PrimeField& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend PrimeField operator-(PrimeField lhs, const PrimeField& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend PrimeField operator*(PrimeField lhs, const PrimeField& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend PrimeField operator/(PrimeField lhs, const PrimeField& rhs) {
        lhs /= rhs;
        return lhs;
    }

    friend constexpr bool operator==(const PrimeField& lhs, const PrimeField& rhs) noexcept =
        default;

    friend constexpr std::strong_ordering operator<=>(
        const PrimeField& lhs, const PrimeField& rhs
    ) noexcept = default;

    static constexpr Integer GetModulus() noexcept {
        return modulus;
    }

private:
    PrimeField ComputeInverse() const;

    constexpr bool IsInvariantSatisfied() const noexcept {
        return value_ >= 0 && value_ < modulus;
    }

    static constexpr Integer Reduce(Integer value) noexcept;

    Integer value_{0};
};

// Implementation details only below here.

template <int64_t modulus>
PrimeField<modulus>::PrimeField(Integer value) : value_(Reduce(value)) {
    assert(IsInvariantSatisfied());
}

template <int64_t modulus>
PrimeField<modulus> PrimeField<modulus>::operator-() const {
    return PrimeField(modulus - value_);
}

template <int64_t modulus>
PrimeField<modulus>& PrimeField<modulus>::operator+=(const PrimeField& addend) {
    value_ += addend.value_;
    if (value_ >= modulus) {
        value_ -= modulus;
    }

    assert(IsInvariantSatisfied());
    return *this;
}

template <int64_t modulus>
PrimeField<modulus>& PrimeField<modulus>::operator-=(const PrimeField& subtrahend) {
    value_ -= subtrahend.value_;
    if (value_ < 0) {
        value_ += modulus;
    }

    assert(IsInvariantSatisfied());
    return *this;
}

template <int64_t modulus>
PrimeField<modulus>& PrimeField<modulus>::operator*=(const PrimeField& factor) {
    value_ = detail::MultiplyModulo(value_, factor.value_, modulus);

    assert(IsInvariantSatisfied());
    return *this;
}

template <int64_t modulus>
PrimeField<modulus>& PrimeField<modulus>::operator/=(const PrimeField& divisor) {
    assert(divisor.value_ != 0 && "The divisor must not be zero");

    *this *= divisor.ComputeInverse();
    assert(IsInvariantSatisfied());
    return *this;
}

template <int64_t modulus>
std::ostream& operator<<(std::ostream& out, const PrimeField<modulus>& number) {
    return out << number.GetValue();
}

template <int64_t modulus>
constexpr typename PrimeField<modulus>::Integer PrimeField<modulus>::Reduce(
    Integer value
) noexcept {
    value %= modulus;

    if (value < 0) {
        value += modulus;
    }

    return value;
}

template <int64_t modulus>
PrimeField<modulus> PrimeField<modulus>::ComputeInverse() const {
    assert(value_ != 0 && "Zero has no multiplicative inverse");

    return PrimeField(detail::ComputeModularExponentiation(value_, modulus - 2, modulus));
}

}  // namespace groebner_basis::fields
