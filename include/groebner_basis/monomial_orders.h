#pragma once

#include <algorithm>

#include "monomial.h"

namespace groebner_basis {
namespace detail {

class TotalDegreeOrder {
public:
    bool operator()(const Monomial& lhs, const Monomial& rhs) const noexcept {
        return lhs.GetTotalDegree() > rhs.GetTotalDegree();
    }
};

class ReverseLexOrder {
public:
    bool operator()(const Monomial& lhs, const Monomial& rhs) const noexcept {
        const auto variable_index_upper_bound =
            std::max(lhs.GetVariableIndexUpperBound(), rhs.GetVariableIndexUpperBound());

        for (auto i = variable_index_upper_bound; i != 0; --i) {
            const auto lhs_degree = lhs.GetDegreeOf(i - 1);
            const auto rhs_degree = rhs.GetDegreeOf(i - 1);

            if (lhs_degree != rhs_degree) {
                return lhs_degree < rhs_degree;
            }
        }

        return false;
    }
};

template <typename PrimaryOrder, typename SecondaryOrder>
class CompositeOrder {
public:
    bool operator()(const Monomial& lhs, const Monomial& rhs) const {
        const PrimaryOrder primary_order{};
        if (primary_order(lhs, rhs)) {
            return true;
        }
        if (primary_order(rhs, lhs)) {
            return false;
        }
        return SecondaryOrder{}(lhs, rhs);
    }
};

}  // namespace detail

class LexicographicOrder {
public:
    bool operator()(const Monomial& lhs, const Monomial& rhs) const noexcept {
        const auto variable_index_upper_bound =
            std::max(lhs.GetVariableIndexUpperBound(), rhs.GetVariableIndexUpperBound());

        for (Monomial::VariableIndex idx = 0; idx != variable_index_upper_bound; ++idx) {
            const auto lhs_degree = lhs.GetDegreeOf(idx);
            const auto rhs_degree = rhs.GetDegreeOf(idx);

            if (lhs_degree != rhs_degree) {
                return lhs_degree > rhs_degree;
            }
        }

        return false;
    }
};

using GradedLexOrder = detail::CompositeOrder<detail::TotalDegreeOrder, LexicographicOrder>;
using GradedReverseLexOrder =
    detail::CompositeOrder<detail::TotalDegreeOrder, detail::ReverseLexOrder>;

}  // namespace groebner_basis
