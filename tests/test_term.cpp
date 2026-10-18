#include <gtest/gtest.h>

#include <sstream>

#include "groebner_basis/fields/rational.h"
#include "groebner_basis/term.h"

namespace groebner_basis {
namespace {

using Coefficient = fields::Rational;
using TermType = Term<Coefficient>;

TEST(TermTest, DefaultConstructionCreatesCanonicalZero) {
    const TermType term;

    EXPECT_TRUE(term.IsZero());
    EXPECT_EQ(term.GetCoefficient(), Coefficient{0});
    EXPECT_TRUE(term.GetMonomial().IsOne());
}

TEST(TermTest, CoefficientAndMonomialConstructorsCreateConstantsAndMonomials) {
    const TermType constant{Coefficient{-3, 2}};
    const TermType monomial{Monomial{1, 0, 2}};

    EXPECT_EQ(constant.GetCoefficient(), Coefficient(-3, 2));
    EXPECT_TRUE(constant.GetMonomial().IsOne());
    EXPECT_EQ(monomial.GetCoefficient(), Coefficient{1});
    EXPECT_EQ(monomial.GetMonomial(), (Monomial{1, 0, 2}));
}

TEST(TermTest, ZeroCoefficientCanonicalizesAwayTheMonomial) {
    const TermType zero{Coefficient{0}, Monomial{4, 3, 2}};

    EXPECT_TRUE(zero.IsZero());
    EXPECT_TRUE(zero.GetMonomial().IsOne());
    EXPECT_EQ(zero, TermType{});
}

TEST(TermTest, DivisibilityUsesMonomialsAndHandlesZeroExplicitly) {
    const TermType divisor{Coefficient{2}, Monomial{1, 0, 2}};
    const TermType dividend{Coefficient{-7, 3}, Monomial{3, 4, 2}};
    const TermType nondivisor{Coefficient{1}, Monomial{4, 0, 2}};
    const TermType zero;

    EXPECT_TRUE(divisor.Divides(dividend));
    EXPECT_TRUE(divisor.Divides(zero));
    EXPECT_FALSE(nondivisor.Divides(dividend));
    EXPECT_FALSE(zero.Divides(dividend));
    EXPECT_FALSE(zero.Divides(zero));
}

TEST(TermTest, UnaryMinusNegatesTheCoefficientAndPreservesTheMonomial) {
    const TermType term{Coefficient{2, 5}, Monomial{1, 3}};

    EXPECT_EQ(-term, (TermType{Coefficient{-2, 5}, Monomial{1, 3}}));
    EXPECT_EQ(-TermType{}, TermType{});
    EXPECT_EQ(term, (TermType{Coefficient{2, 5}, Monomial{1, 3}}));
}

TEST(TermTest, CoefficientMultiplicationHandlesIdentitySignAndZero) {
    const TermType term{Coefficient{2, 3}, Monomial{1, 2}};

    EXPECT_EQ(term * Coefficient{1}, term);
    EXPECT_EQ(Coefficient{-3} * term, (TermType{Coefficient{-2}, Monomial{1, 2}}));
    EXPECT_EQ(term * Coefficient{0}, TermType{});
    EXPECT_EQ(term, (TermType{Coefficient{2, 3}, Monomial{1, 2}}));
}

TEST(TermTest, MonomialMultiplicationAddsVariableDegrees) {
    const TermType term{Coefficient{-2}, Monomial{1, 0, 2}};
    const Monomial factor{0, 3, 1, 4};

    EXPECT_EQ(term * factor, (TermType{Coefficient{-2}, Monomial{1, 3, 3, 4}}));
    EXPECT_EQ(factor * term, (TermType{Coefficient{-2}, Monomial{1, 3, 3, 4}}));
    EXPECT_EQ(term * Monomial{}, term);
}

TEST(TermTest, TermMultiplicationCombinesCoefficientsAndMonomials) {
    const TermType lhs{Coefficient{2, 3}, Monomial{1, 0, 2}};
    const TermType rhs{Coefficient{-9, 4}, Monomial{0, 3, 1}};

    EXPECT_EQ(lhs * rhs, (TermType{Coefficient{-3, 2}, Monomial{1, 3, 3}}));
    EXPECT_EQ(lhs * TermType{}, TermType{});
    EXPECT_EQ(TermType{} * rhs, TermType{});
    EXPECT_EQ(lhs, (TermType{Coefficient{2, 3}, Monomial{1, 0, 2}}));
    EXPECT_EQ(rhs, (TermType{Coefficient{-9, 4}, Monomial{0, 3, 1}}));
}

TEST(TermTest, MultiplicationAssignmentsReturnTheLeftOperand) {
    TermType term{Coefficient{2}, Monomial{1}};

    EXPECT_EQ(&(term *= Coefficient{-1}), &term);
    EXPECT_EQ(&(term *= Monomial{0, 2}), &term);
    const TermType factor{Coefficient{3}, Monomial{1}};
    EXPECT_EQ(&(term *= factor), &term);
    EXPECT_EQ(term, (TermType{Coefficient{-6}, Monomial{2, 2}}));
}

TEST(TermTest, DivisionByCoefficientMonomialAndTermProducesExpectedQuotients) {
    const TermType dividend{Coefficient{-6, 5}, Monomial{3, 4, 2}};

    EXPECT_EQ((dividend / Coefficient{3, 5}), (TermType{Coefficient{-2}, Monomial{3, 4, 2}}));
    EXPECT_EQ((dividend / Monomial{1, 4, 1}), (TermType{Coefficient{-6, 5}, Monomial{2, 0, 1}}));
    EXPECT_EQ(
        (dividend / TermType{Coefficient{3, 2}, Monomial{1, 4, 1}}),
        (TermType{Coefficient{-4, 5}, Monomial{2, 0, 1}})
    );
    EXPECT_EQ(TermType{} / Coefficient{2}, TermType{});
    EXPECT_EQ((TermType{} / Monomial{9, 9}), TermType{});
}

TEST(TermTest, DivisionAssignmentsReturnTheLeftOperand) {
    TermType term{Coefficient{6}, Monomial{3, 2}};

    EXPECT_EQ(&(term /= Coefficient{2}), &term);
    EXPECT_EQ(&(term /= Monomial{1, 1}), &term);
    const TermType divisor{Coefficient{3}, Monomial{1}};
    EXPECT_EQ(&(term /= divisor), &term);
    EXPECT_EQ(term, (TermType{Coefficient{1}, Monomial{1, 1}}));
}

#ifndef NDEBUG
TEST(TermDeathTest, DivisionRejectsInvalidDivisors) {
    EXPECT_DEATH(
        {
            TermType dividend(Coefficient(1), Monomial{1});
            dividend /= Coefficient{0};
        },
        "nonzero coefficient"
    );
    EXPECT_DEATH(
        {
            TermType dividend(Coefficient(1), Monomial{1, 1});
            dividend /= Monomial{2};
        },
        "divisor monomial divides dividend"
    );
    EXPECT_DEATH(
        {
            TermType dividend(Coefficient(1), Monomial{1});
            dividend /= TermType{};
        },
        "divisor divides dividend"
    );
}
#endif

TEST(TermTest, EqualityComparesBothCoefficientAndMonomial) {
    const TermType term{Coefficient{2, 3}, Monomial{1, 2}};

    EXPECT_EQ(term, (TermType{Coefficient{4, 6}, Monomial{1, 2}}));
    EXPECT_NE(term, (TermType{Coefficient{2, 3}, Monomial{2, 1}}));
    EXPECT_NE(term, (TermType{Coefficient{3, 2}, Monomial{1, 2}}));
    EXPECT_EQ((TermType{Coefficient{0}, Monomial{7}}), TermType{});
}

TEST(TermTest, StreamOutputHandlesEveryCoefficientSpecialCase) {
    std::ostringstream output;

    output << TermType{} << " | " << TermType{Coefficient{-3, 2}} << " | "
           << TermType{Monomial{1, 0, 2}} << " | " << TermType{Coefficient{-1}, Monomial{0, 1}}
           << " | " << TermType{Coefficient{2, 3}, Monomial{1, 2}};

    EXPECT_EQ(output.str(), "0 | -3/2 | x_0*x_2^2 | -x_1 | 2/3*x_0*x_1^2");
}

}  // namespace
}  // namespace groebner_basis
