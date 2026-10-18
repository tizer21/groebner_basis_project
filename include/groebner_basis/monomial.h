#pragma once

#include <stddef.h>

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <numeric>
#include <ostream>
#include <utility>
#include <vector>

namespace groebner_basis {

class Monomial {
public:
    using Degree = unsigned int;
    using VariableIndex = size_t;

    Monomial() = default;
    explicit Monomial(std::initializer_list<Degree> degrees);
    explicit Monomial(std::vector<Degree> degrees);

    Monomial(const Monomial&) = default;
    Monomial(Monomial&&) noexcept;

    Monomial& operator=(const Monomial&) = default;
    Monomial& operator=(Monomial&&) noexcept;

    bool Divides(const Monomial& dividend) const noexcept;

    Degree GetTotalDegree() const noexcept {
        return total_degree_;
    }

    Degree GetDegreeOf(VariableIndex idx) const noexcept {
        return (idx < degrees_.size()) ? degrees_[idx] : Degree{0};
    }

    VariableIndex GetVariableIndexUpperBound() const noexcept {
        return degrees_.size();
    }

    bool IsOne() const noexcept {
        return degrees_.empty();
    }

    Monomial& operator*=(const Monomial& factor);
    Monomial& operator/=(const Monomial& divisor);

    friend Monomial operator*(Monomial lhs, const Monomial& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Monomial operator/(Monomial lhs, const Monomial& rhs) {
        lhs /= rhs;
        return lhs;
    }

    friend bool operator==(const Monomial& lhs, const Monomial& rhs) noexcept {
        return lhs.degrees_ == rhs.degrees_;
    }

    friend void swap(Monomial& lhs, Monomial& rhs) noexcept {
        lhs.Swap(rhs);
    }

private:
    using DegreeContainer = std::vector<Degree>;

    static DegreeContainer TrimTrailingZeros(DegreeContainer&& degrees) noexcept {
        while (!degrees.empty() && degrees.back() == Degree{0}) {
            degrees.pop_back();
        }
        return degrees;
    }

    static Degree CalculateTotalDegree(const DegreeContainer& degrees) noexcept {
        return std::accumulate(degrees.begin(), degrees.end(), Degree{0});
    }

    bool IsInvariantSatisfied() const noexcept {
        const bool is_canonical = degrees_.empty() || degrees_.back() != Degree{0};
        return is_canonical && total_degree_ == CalculateTotalDegree(degrees_);
    }

    void Swap(Monomial& other) noexcept {
        degrees_.swap(other.degrees_);
        std::swap(total_degree_, other.total_degree_);
    }

    DegreeContainer degrees_;
    Degree total_degree_{0};
};

inline std::ostream& operator<<(std::ostream& os, const Monomial& monomial);

inline Monomial ComputeLeastCommonMultiple(const Monomial& lhs, const Monomial& rhs);

// Implementation details only below here.

inline Monomial::Monomial(std::initializer_list<Degree> degrees)
    : Monomial(DegreeContainer(degrees)) {}

inline Monomial::Monomial(std::vector<Degree> degrees)
    : degrees_(TrimTrailingZeros(std::move(degrees))),
      total_degree_(CalculateTotalDegree(degrees_)) {
    assert(IsInvariantSatisfied());
}

inline Monomial::Monomial(Monomial&& other) noexcept : Monomial() {
    Swap(other);

    assert(IsInvariantSatisfied());
    assert(other.IsInvariantSatisfied());
}

inline Monomial& Monomial::operator=(Monomial&& other) noexcept {
    Monomial temporary(std::move(other));

    Swap(temporary);

    assert(IsInvariantSatisfied());
    assert(other.IsInvariantSatisfied());
    return *this;
}

inline bool Monomial::Divides(const Monomial& dividend) const noexcept {
    if (degrees_.size() > dividend.degrees_.size()) {
        return false;
    }
    for (VariableIndex idx = 0; idx != degrees_.size(); ++idx) {
        if (degrees_[idx] > dividend.degrees_[idx]) {
            return false;
        }
    }
    return true;
}

inline Monomial& Monomial::operator*=(const Monomial& factor) {
    const auto factor_size = factor.degrees_.size();
    if (degrees_.size() < factor_size) {
        degrees_.resize(factor_size);
    }
    for (VariableIndex idx = 0; idx != factor_size; ++idx) {
        degrees_[idx] += factor.degrees_[idx];
    }

    total_degree_ += factor.total_degree_;

    assert(IsInvariantSatisfied());
    return *this;
}

inline Monomial& Monomial::operator/=(const Monomial& divisor) {
    assert(divisor.Divides(*this) && "Monomial division requires that divisor divides dividend");

    for (VariableIndex idx = 0; idx != divisor.degrees_.size(); ++idx) {
        degrees_[idx] -= divisor.degrees_[idx];
    }
    degrees_ = TrimTrailingZeros(std::move(degrees_));

    total_degree_ -= divisor.total_degree_;

    assert(IsInvariantSatisfied());
    return *this;
}

inline std::ostream& operator<<(std::ostream& os, const Monomial& monomial) {
    if (monomial.IsOne()) {
        return os << '1';
    }

    bool is_first = true;
    for (Monomial::VariableIndex idx = 0; idx < monomial.GetVariableIndexUpperBound(); ++idx) {
        const auto degree = monomial.GetDegreeOf(idx);
        if (degree == Monomial::Degree{0}) {
            continue;
        }

        if (!is_first) {
            os << '*';
        }
        os << "x_" << idx;
        if (degree > Monomial::Degree{1}) {
            os << '^' << degree;
        }

        is_first = false;
    }

    return os;
}

inline Monomial ComputeLeastCommonMultiple(const Monomial& lhs, const Monomial& rhs) {
    const auto variable_index_upper_bound =
        std::max(lhs.GetVariableIndexUpperBound(), rhs.GetVariableIndexUpperBound());
    std::vector<Monomial::Degree> degrees(variable_index_upper_bound);
    for (Monomial::VariableIndex idx = 0; idx != variable_index_upper_bound; ++idx) {
        degrees[idx] = std::max(lhs.GetDegreeOf(idx), rhs.GetDegreeOf(idx));
    }

    return Monomial{std::move(degrees)};
}

}  // namespace groebner_basis
