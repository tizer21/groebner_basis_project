#pragma once

#include <cassert>
#include <initializer_list>
#include <utility>
#include <vector>

#include "algorithms.h"
#include "polynomial.h"

namespace groebner_basis {

template <typename CoefficientType, typename MonomialOrder>
class PolynomialIdeal {
public:
    using PolynomialType = Polynomial<CoefficientType, MonomialOrder>;
    using PolynomialContainer = std::vector<PolynomialType>;

    PolynomialIdeal() = default;

    explicit PolynomialIdeal(std::initializer_list<PolynomialType> generators)
        : PolynomialIdeal(PolynomialContainer(generators)) {}

    explicit PolynomialIdeal(std::vector<PolynomialType> generators)
        : reduced_basis_(algorithms::ComputeReducedGroebnerBasis(std::move(generators))) {}

    const PolynomialContainer& GetReducedGroebnerBasis() const noexcept {
        return reduced_basis_;
    }

    bool IsZero() const noexcept {
        return reduced_basis_.empty();
    }

    bool IsUnit() const {
        return reduced_basis_.size() == 1 && reduced_basis_.front().IsOne();
    }

    bool ContainsPolynomial(const PolynomialType& polynomial) const {
        return ComputeNormalForm(polynomial).IsZero();
    }

    PolynomialType ComputeNormalForm(const PolynomialType& polynomial) const {
        return algorithms::ReducePolynomial(polynomial, reduced_basis_);
    }

    PolynomialIdeal& operator+=(const PolynomialIdeal& other) {
        if (this == &other) {
            return *this;
        }

        auto basis_generators = reduced_basis_;
        basis_generators.insert(
            basis_generators.end(), other.reduced_basis_.begin(), other.reduced_basis_.end()
        );

        auto new_reduced_basis =
            algorithms::ComputeReducedGroebnerBasis(std::move(basis_generators));
        reduced_basis_ = std::move(new_reduced_basis);

        return *this;
    }

    PolynomialIdeal& operator*=(const PolynomialIdeal& other) {
        PolynomialContainer basis_generators;
        basis_generators.reserve(reduced_basis_.size() * other.reduced_basis_.size());
        for (const auto& lhs : reduced_basis_) {
            for (const auto& rhs : other.reduced_basis_) {
                basis_generators.push_back(lhs * rhs);
            }
        }

        auto new_reduced_basis =
            algorithms::ComputeReducedGroebnerBasis(std::move(basis_generators));
        reduced_basis_ = std::move(new_reduced_basis);

        return *this;
    }

    bool operator==(const PolynomialIdeal& other) const {
        return reduced_basis_ == other.reduced_basis_;
    }

    friend PolynomialIdeal operator+(PolynomialIdeal lhs, const PolynomialIdeal& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend PolynomialIdeal operator*(PolynomialIdeal lhs, const PolynomialIdeal& rhs) {
        lhs *= rhs;
        return lhs;
    }

private:
    PolynomialContainer reduced_basis_;
};

}  // namespace groebner_basis
