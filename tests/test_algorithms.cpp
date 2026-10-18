#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "groebner_basis/algorithms.h"
#include "groebner_basis/fields/prime_field.h"
#include "groebner_basis/fields/rational.h"
#include "groebner_basis/monomial_orders.h"

namespace groebner_basis::algorithms {
namespace {

using Coefficient = fields::Rational;
using TermType = Term<Coefficient>;
using PolynomialType = Polynomial<Coefficient, LexicographicOrder>;
using Basis = PolynomialBasis<Coefficient, LexicographicOrder>;

template <typename Field, typename Order>
void ExpectGroebnerBasis(const PolynomialBasis<Field, Order>& basis) {
    for (const auto& polynomial : basis) {
        EXPECT_FALSE(polynomial.IsZero());
    }

    for (std::size_t lhs = 0; lhs < basis.size(); ++lhs) {
        for (std::size_t rhs = lhs + 1; rhs < basis.size(); ++rhs) {
            SCOPED_TRACE("S-pair (" + std::to_string(lhs) + ", " + std::to_string(rhs) + ")");
            const auto s_polynomial = ComputeSPolynomial(basis[lhs], basis[rhs]);
            EXPECT_TRUE(ReducePolynomial(s_polynomial, basis).IsZero());
        }
    }
}

template <typename Field, typename Order>
void ExpectPolynomialDivision(
    const Polynomial<Field, Order>& dividend, const PolynomialBasis<Field, Order>& divisors,
    const PolynomialBasis<Field, Order>& expected_quotients,
    const Polynomial<Field, Order>& expected_remainder
) {
    const auto result = DividePolynomial(dividend, divisors);

    ASSERT_EQ(result.quotients.size(), divisors.size());
    EXPECT_EQ(result.quotients, expected_quotients);
    EXPECT_EQ(result.remainder, expected_remainder);
    EXPECT_EQ(result.remainder, ReducePolynomial(dividend, divisors));

    auto reconstructed_dividend = result.remainder;
    for (std::size_t idx = 0; idx != divisors.size(); ++idx) {
        reconstructed_dividend += result.quotients[idx] * divisors[idx];
    }
    EXPECT_EQ(reconstructed_dividend, dividend);
}

TEST(BuchbergerDetailTest, PairConstructionCanonicalizesIndices) {
    EXPECT_EQ(detail::MakeBuchbergerPair(1, 4), std::make_pair(1UL, 4UL));
    EXPECT_EQ(detail::MakeBuchbergerPair(4, 1), std::make_pair(1UL, 4UL));
}

#ifndef NDEBUG
TEST(BuchbergerDetailDeathTest, PairConstructionRejectsEqualIndices) {
    EXPECT_DEATH((void)detail::MakeBuchbergerPair(2, 2), "");
}
#endif

TEST(BuchbergerDetailTest, PairSetContainsEveryUnorderedPairExactlyOnce) {
    const Basis basis{
        PolynomialType{TermType{Coefficient{1}, Monomial{3}}},
        PolynomialType{TermType{Coefficient{1}, Monomial{2}}},
        PolynomialType{TermType{Coefficient{1}, Monomial{1}}},
        PolynomialType{TermType{Coefficient{1}}},
    };
    const detail::BuchbergerPairSet expected{
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3},
    };

    EXPECT_EQ(detail::BuildBuchbergerPairs(basis), expected);
    EXPECT_TRUE(detail::BuildBuchbergerPairs(Basis{}).empty());
}

TEST(BuchbergerDetailTest, ProductCriterionDetectsRelativelyPrimeLeadingMonomials) {
    const PolynomialType x_squared{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{1}},
    };
    const PolynomialType y_squared{
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{1}},
    };
    const PolynomialType x_y{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{1}},
    };

    EXPECT_TRUE(detail::MeetsBuchbergerProductCriterion(x_squared, y_squared));
    EXPECT_FALSE(detail::MeetsBuchbergerProductCriterion(x_squared, x_y));
}

TEST(BuchbergerDetailTest, ChainCriterionRequiresBothIntermediatePairsToBeProcessed) {
    const Basis basis{
        PolynomialType{TermType{Coefficient{1}, Monomial{2}}},
        PolynomialType{TermType{Coefficient{1}, Monomial{0, 2}}},
        PolynomialType{TermType{Coefficient{1}, Monomial{1, 1}}},
    };

    EXPECT_TRUE(detail::MeetsBuchbergerChainCriterion(0, 1, basis, detail::BuchbergerPairSet{}));

    const detail::BuchbergerPairSet pending{{0, 2}};
    EXPECT_FALSE(detail::MeetsBuchbergerChainCriterion(0, 1, basis, pending));
}

TEST(DividePolynomialTest, ZeroDividendAndEmptyDivisorListAreHandled) {
    const PolynomialType zero;
    const PolynomialType polynomial{
        TermType{Coefficient{2}, Monomial{2}},
        TermType{Coefficient{-3}, Monomial{0, 1}},
        TermType{Coefficient{1}},
    };

    ExpectPolynomialDivision(zero, Basis{}, Basis{}, zero);
    ExpectPolynomialDivision(polynomial, Basis{}, Basis{}, polynomial);
    ExpectPolynomialDivision(zero, Basis{polynomial}, Basis{PolynomialType{}}, zero);
}

TEST(DividePolynomialTest, ZeroDivisorsAreSkippedAndRetainTheirQuotientPositions) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType divisor{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType quotient{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };

    ExpectPolynomialDivision(
        dividend, Basis{PolynomialType{}, divisor}, Basis{PolynomialType{}, quotient},
        PolynomialType{}
    );
}

TEST(DividePolynomialTest, DivisionUsesCoefficientRatios) {
    const PolynomialType dividend{
        TermType{Coefficient{2}, Monomial{2}},
        TermType{Coefficient{3}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType divisor{
        TermType{Coefficient{4}, Monomial{1}},
        TermType{Coefficient{2}},
    };
    const PolynomialType quotient{
        TermType{Coefficient{1, 2}, Monomial{1}},
        TermType{Coefficient{1, 2}},
    };

    ExpectPolynomialDivision(dividend, Basis{divisor}, Basis{quotient}, PolynomialType{});
}

TEST(DividePolynomialTest, CompleteUnivariateDivisionProducesQuotientAndRemainder) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{3}},
        TermType{Coefficient{2}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType divisor{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType expected_quotient{
        TermType{Coefficient{1, 2}, Monomial{2}},
        TermType{Coefficient{3, 4}, Monomial{1}},
        TermType{Coefficient{1, 8}},
    };
    const PolynomialType expected_remainder{TermType{Coefficient{7, 8}}};

    ExpectPolynomialDivision(
        dividend, Basis{divisor}, Basis{expected_quotient}, expected_remainder
    );
}

TEST(DividePolynomialTest, TextbookDivisionProducesBothQuotientsAndConstantRemainder) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{1, 2}},
        TermType{Coefficient{1}},
    };
    const PolynomialType x_y_plus_one{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType y_plus_one{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType y{TermType{Coefficient{1}, Monomial{0, 1}}};
    const PolynomialType negative_one{TermType{Coefficient{-1}}};
    const PolynomialType expected_remainder{TermType{Coefficient{2}}};

    ExpectPolynomialDivision(
        dividend, Basis{x_y_plus_one, y_plus_one}, Basis{y, negative_one}, expected_remainder
    );
}

TEST(DividePolynomialTest, TextbookMultivariateDivisionProducesBothQuotients) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{2, 1}},
        TermType{Coefficient{1}, Monomial{1, 2}},
        TermType{Coefficient{1}, Monomial{0, 2}},
    };
    const PolynomialType x_y_minus_one{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType y_squared_minus_one{
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType x_plus_y{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };
    const PolynomialType one{TermType{Coefficient{1}}};
    const PolynomialType expected_remainder{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{1}},
    };

    ExpectPolynomialDivision(
        dividend, Basis{x_y_minus_one, y_squared_minus_one}, Basis{x_plus_y, one},
        expected_remainder
    );
}

TEST(DividePolynomialTest, ReversedTextbookDivisorsProduceDifferentQuotients) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{2, 1}},
        TermType{Coefficient{1}, Monomial{1, 2}},
        TermType{Coefficient{1}, Monomial{0, 2}},
    };
    const PolynomialType y_squared_minus_one{
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType x_y_minus_one{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType x_plus_one{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType x{TermType{Coefficient{1}, Monomial{1}}};
    const PolynomialType expected_remainder{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{1}},
    };

    ExpectPolynomialDivision(
        dividend, Basis{y_squared_minus_one, x_y_minus_one}, Basis{x_plus_one, x},
        expected_remainder
    );
}

TEST(DividePolynomialTest, TextbookLexDivisionLeavesUnusedDivisorQuotientAtZero) {
    const PolynomialType dividend{TermType{Coefficient{1}, Monomial{5, 1}}};
    const PolynomialType x_squared_y_minus_y_squared{
        TermType{Coefficient{1}, Monomial{2, 1}},
        TermType{Coefficient{-1}, Monomial{0, 2}},
    };
    const PolynomialType x_fourth_y_squared_minus_y_squared{
        TermType{Coefficient{1}, Monomial{4, 2}},
        TermType{Coefficient{-1}, Monomial{0, 2}},
    };
    const PolynomialType expected_quotient{
        TermType{Coefficient{1}, Monomial{3}},
        TermType{Coefficient{1}, Monomial{1, 1}},
    };
    const PolynomialType expected_remainder{
        TermType{Coefficient{1}, Monomial{1, 3}},
    };

    ExpectPolynomialDivision(
        dividend, Basis{x_squared_y_minus_y_squared, x_fourth_y_squared_minus_y_squared},
        Basis{expected_quotient, PolynomialType{}}, expected_remainder
    );
}

TEST(DividePolynomialTest, MultivariateQuotientsDependOnDivisorOrder) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{1, 2}},
        TermType{Coefficient{-1}, Monomial{1}},
    };
    const PolynomialType x_y_minus_one{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType y_squared_minus_one{
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType x{TermType{Coefficient{1}, Monomial{1}}};
    const PolynomialType y{TermType{Coefficient{1}, Monomial{0, 1}}};
    const PolynomialType expected_first_remainder{
        TermType{Coefficient{-1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };

    ExpectPolynomialDivision(
        dividend, Basis{x_y_minus_one, y_squared_minus_one}, Basis{y, PolynomialType{}},
        expected_first_remainder
    );
    ExpectPolynomialDivision(
        dividend, Basis{y_squared_minus_one, x_y_minus_one}, Basis{x, PolynomialType{}},
        PolynomialType{}
    );
}

TEST(ReducePolynomialTest, ZeroDividendAndEmptyDivisorListAreHandled) {
    const PolynomialType zero;
    const PolynomialType polynomial{
        TermType{Coefficient{2}, Monomial{2}},
        TermType{Coefficient{-3}, Monomial{0, 1}},
        TermType{Coefficient{1}},
    };

    EXPECT_TRUE(ReducePolynomial(zero, Basis{}).IsZero());
    EXPECT_EQ(ReducePolynomial(polynomial, Basis{}), polynomial);
}

TEST(ReducePolynomialTest, ZeroDivisorsAreSkipped) {
    const PolynomialType polynomial{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType divisor{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{-1}},
    };

    EXPECT_TRUE(ReducePolynomial(polynomial, Basis{PolynomialType{}, divisor}).IsZero());
}

TEST(ReducePolynomialTest, DivisionUsesCoefficientRatios) {
    const PolynomialType dividend{
        TermType{Coefficient{2}, Monomial{2}},
        TermType{Coefficient{3}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType divisor{
        TermType{Coefficient{4}, Monomial{1}},
        TermType{Coefficient{2}},
    };

    EXPECT_TRUE(ReducePolynomial(dividend, Basis{divisor}).IsZero());
}

TEST(ReducePolynomialTest, OneLeadingMonomialReductionProducesTheExpectedPolynomial) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{3}},
        TermType{Coefficient{2}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType divisor{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType expected{
        TermType{Coefficient{3, 2}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };

    const auto leading_term_factor = dividend.GetLeadingTerm() / divisor.GetLeadingTerm();
    const auto reduced_once = dividend - leading_term_factor * divisor;

    EXPECT_EQ(reduced_once, expected);
}

TEST(ReducePolynomialTest, CompleteUnivariateReductionContinuesToConstantRemainder) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{3}},
        TermType{Coefficient{2}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType divisor{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType expected{TermType{Coefficient{7, 8}}};

    EXPECT_EQ(ReducePolynomial(dividend, Basis{divisor}), expected);
}

TEST(ReducePolynomialTest, TextbookDivisionProducesConstantRemainderAndIdentity) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{1, 2}},
        TermType{Coefficient{1}},
    };
    const PolynomialType x_y_plus_one{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType y_plus_one{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType y{TermType{Coefficient{1}, Monomial{0, 1}}};
    const PolynomialType expected_remainder{TermType{Coefficient{2}}};
    const Basis divisors{x_y_plus_one, y_plus_one};

    EXPECT_EQ(dividend, y * x_y_plus_one - y_plus_one + expected_remainder);
    EXPECT_EQ(ReducePolynomial(dividend, divisors), expected_remainder);
    EXPECT_EQ(ReducePolynomial(expected_remainder, divisors), expected_remainder);
}

TEST(ReducePolynomialTest, TextbookMultivariateDivisionProducesExpectedRemainderAndIdentity) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{2, 1}},
        TermType{Coefficient{1}, Monomial{1, 2}},
        TermType{Coefficient{1}, Monomial{0, 2}},
    };
    const PolynomialType x_y_minus_one{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType y_squared_minus_one{
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType expected{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType x_plus_y{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };
    const Basis divisors{x_y_minus_one, y_squared_minus_one};

    EXPECT_EQ(dividend, x_plus_y * x_y_minus_one + y_squared_minus_one + expected);
    EXPECT_EQ(ReducePolynomial(dividend, divisors), expected);
    EXPECT_EQ(ReducePolynomial(expected, divisors), expected);
}

TEST(ReducePolynomialTest, ReversedTextbookDivisorsProduceDifferentRemainderAndIdentity) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{2, 1}},
        TermType{Coefficient{1}, Monomial{1, 2}},
        TermType{Coefficient{1}, Monomial{0, 2}},
    };
    const PolynomialType y_squared_minus_one{
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType x_y_minus_one{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType x_plus_one{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType x{TermType{Coefficient{1}, Monomial{1}}};
    const PolynomialType expected_remainder{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const Basis divisors{y_squared_minus_one, x_y_minus_one};

    EXPECT_EQ(dividend, x_plus_one * y_squared_minus_one + x * x_y_minus_one + expected_remainder);
    EXPECT_EQ(ReducePolynomial(dividend, divisors), expected_remainder);
    EXPECT_EQ(ReducePolynomial(expected_remainder, divisors), expected_remainder);
}

TEST(ReducePolynomialTest, TextbookLexDivisionProducesXYCubedRemainder) {
    const PolynomialType dividend{TermType{Coefficient{1}, Monomial{5, 1}}};
    const PolynomialType x_squared_y_minus_y_squared{
        TermType{Coefficient{1}, Monomial{2, 1}},
        TermType{Coefficient{-1}, Monomial{0, 2}},
    };
    const PolynomialType x_fourth_y_squared_minus_y_squared{
        TermType{Coefficient{1}, Monomial{4, 2}},
        TermType{Coefficient{-1}, Monomial{0, 2}},
    };
    const PolynomialType quotient{
        TermType{Coefficient{1}, Monomial{3}},
        TermType{Coefficient{1}, Monomial{1, 1}},
    };
    const PolynomialType expected_remainder{TermType{Coefficient{1}, Monomial{1, 3}}};
    const Basis divisors{
        x_squared_y_minus_y_squared,
        x_fourth_y_squared_minus_y_squared,
    };

    EXPECT_EQ(dividend, quotient * x_squared_y_minus_y_squared + expected_remainder);
    EXPECT_EQ(ReducePolynomial(dividend, divisors), expected_remainder);
    EXPECT_EQ(ReducePolynomial(expected_remainder, divisors), expected_remainder);
}

TEST(ReducePolynomialTest, MultivariateRemainderCanDependOnDivisorOrder) {
    const PolynomialType dividend{
        TermType{Coefficient{1}, Monomial{1, 2}},
        TermType{Coefficient{-1}, Monomial{1}},
    };
    const PolynomialType x_y_minus_one{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType y_squared_minus_one{
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType expected_first_order{
        TermType{Coefficient{-1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };
    const PolynomialType x{TermType{Coefficient{1}, Monomial{1}}};
    const PolynomialType y{TermType{Coefficient{1}, Monomial{0, 1}}};
    const Basis first_order{x_y_minus_one, y_squared_minus_one};
    const Basis reversed_order{y_squared_minus_one, x_y_minus_one};

    EXPECT_EQ(dividend, y * x_y_minus_one + expected_first_order);
    EXPECT_EQ(ReducePolynomial(dividend, first_order), expected_first_order);
    EXPECT_EQ(ReducePolynomial(expected_first_order, first_order), expected_first_order);

    EXPECT_EQ(dividend, x * y_squared_minus_one);
    EXPECT_TRUE(ReducePolynomial(dividend, reversed_order).IsZero());
}

TEST(ComputeSPolynomialTest, LeadingTermsCancelIncludingTheirCoefficients) {
    const PolynomialType lhs{
        TermType{Coefficient{2}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };
    const PolynomialType rhs{
        TermType{Coefficient{3}, Monomial{1, 1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType expected{
        TermType{Coefficient{-1, 3}, Monomial{1}},
        TermType{Coefficient{1, 2}, Monomial{0, 2}},
    };

    EXPECT_EQ(ComputeSPolynomial(lhs, rhs), expected);
    EXPECT_EQ(ComputeSPolynomial(rhs, lhs), -expected);
    EXPECT_TRUE(ComputeSPolynomial(lhs, lhs).IsZero());
}

TEST(ComputeSPolynomialTest, TextbookGradedLexExampleUsesLcmAndLeadingCoefficients) {
    using GradedPolynomial = Polynomial<Coefficient, GradedLexOrder>;

    const GradedPolynomial f{
        TermType{Coefficient{1}, Monomial{3, 2}},
        TermType{Coefficient{-1}, Monomial{2, 3}},
        TermType{Coefficient{1}, Monomial{1}},
    };
    const GradedPolynomial g{
        TermType{Coefficient{3}, Monomial{4, 1}},
        TermType{Coefficient{1}, Monomial{0, 2}},
    };
    const GradedPolynomial expected{
        TermType{Coefficient{-1}, Monomial{3, 3}},
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-1, 3}, Monomial{0, 3}},
    };
    const TermType x_factor{Coefficient{1}, Monomial{1}};
    const TermType one_third_y_factor{
        Coefficient{1, 3},
        Monomial{0, 1},
    };

    EXPECT_EQ(f.GetLeadingMonomial(), (Monomial{3, 2}));
    EXPECT_EQ(g.GetLeadingMonomial(), (Monomial{4, 1}));
    EXPECT_EQ(
        ComputeLeastCommonMultiple(f.GetLeadingMonomial(), g.GetLeadingMonomial()), (Monomial{4, 2})
    );

    EXPECT_EQ(ComputeSPolynomial(f, g), x_factor * f - one_third_y_factor * g);
    EXPECT_EQ(ComputeSPolynomial(f, g), expected);
}

TEST(ComputeSPolynomialTest, TextbookLexExampleUsesVariablePriorityYThenZThenX) {
    // Degree coordinates are (y, z, x), so lexicographic order is y > z > x.
    const PolynomialType y_minus_x_squared{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
    };
    const PolynomialType z_minus_x_cubed{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 3}},
    };
    const PolynomialType expected{
        TermType{Coefficient{-1}, Monomial{0, 1, 2}},
        TermType{Coefficient{1}, Monomial{1, 0, 3}},
    };
    const TermType z_factor{Coefficient{1}, Monomial{0, 1}};
    const TermType y_factor{Coefficient{1}, Monomial{1}};
    const TermType x_cubed_factor{Coefficient{1}, Monomial{0, 0, 3}};
    const TermType negative_x_squared_factor{
        Coefficient{-1},
        Monomial{0, 0, 2},
    };
    const Basis generators{y_minus_x_squared, z_minus_x_cubed};

    EXPECT_EQ(y_minus_x_squared.GetLeadingMonomial(), (Monomial{1}));
    EXPECT_EQ(z_minus_x_cubed.GetLeadingMonomial(), (Monomial{0, 1}));
    EXPECT_EQ(
        ComputeLeastCommonMultiple(
            y_minus_x_squared.GetLeadingMonomial(), z_minus_x_cubed.GetLeadingMonomial()
        ),
        (Monomial{1, 1})
    );

    EXPECT_EQ(
        ComputeSPolynomial(y_minus_x_squared, z_minus_x_cubed),
        z_factor * y_minus_x_squared - y_factor * z_minus_x_cubed
    );
    EXPECT_EQ(ComputeSPolynomial(y_minus_x_squared, z_minus_x_cubed), expected);
    EXPECT_EQ(
        expected, x_cubed_factor * y_minus_x_squared + negative_x_squared_factor * z_minus_x_cubed
    );
    EXPECT_TRUE(ReducePolynomial(expected, generators).IsZero());
    EXPECT_EQ(ComputeGroebnerBasis(generators), generators);
    EXPECT_EQ(ComputeReducedGroebnerBasis(generators), generators);
    ExpectGroebnerBasis(generators);
}

TEST(BuchbergerAlgorithmTest, TwistedCubicGeneratorsAreNotAGroebnerBasisWhenXIsHighest) {
    // Degree coordinates are (x, y, z), so lexicographic order is x > y > z.
    const PolynomialType y_minus_x_squared{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{2}},
    };
    const PolynomialType z_minus_x_cubed{
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{-1}, Monomial{3}},
    };
    const PolynomialType expected_s_polynomial{
        TermType{Coefficient{-1}, Monomial{1, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
    };
    const Basis generators{y_minus_x_squared, z_minus_x_cubed};

    EXPECT_EQ(y_minus_x_squared.GetLeadingMonomial(), (Monomial{2}));
    EXPECT_EQ(z_minus_x_cubed.GetLeadingMonomial(), (Monomial{3}));

    const auto s_polynomial = ComputeSPolynomial(y_minus_x_squared, z_minus_x_cubed);
    EXPECT_EQ(s_polynomial, expected_s_polynomial);
    EXPECT_EQ(ReducePolynomial(s_polynomial, generators), s_polynomial);
    EXPECT_FALSE(y_minus_x_squared.GetLeadingMonomial().Divides(s_polynomial.GetLeadingMonomial()));
    EXPECT_FALSE(z_minus_x_cubed.GetLeadingMonomial().Divides(s_polynomial.GetLeadingMonomial()));

    const auto completed_basis = ComputeGroebnerBasis(generators);
    EXPECT_GT(completed_basis.size(), generators.size());
    ExpectGroebnerBasis(completed_basis);
}

TEST(BuchbergerAlgorithmTest, ProductCriterionDoesNotRequireTheFixedDivisionRemainderToBeZero) {
    using GradedPolynomial = Polynomial<Coefficient, GradedLexOrder>;
    using GradedBasis = PolynomialBasis<Coefficient, GradedLexOrder>;

    const GradedPolynomial y_z_plus_y{
        TermType{Coefficient{1}, Monomial{0, 1, 1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };
    const GradedPolynomial x_cubed_plus_y{
        TermType{Coefficient{1}, Monomial{3}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };
    const GradedPolynomial z_fourth{TermType{Coefficient{1}, Monomial{0, 0, 4}}};
    const GradedPolynomial expected_s_polynomial{TermType{Coefficient{1}, Monomial{0, 1, 4}}};
    const GradedPolynomial quotient{
        TermType{Coefficient{1}, Monomial{0, 0, 3}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{-1}},
    };
    const GradedPolynomial expected_fixed_remainder{TermType{Coefficient{1}, Monomial{0, 1}}};
    const GradedBasis ordered_divisors{
        y_z_plus_y,
        x_cubed_plus_y,
        z_fourth,
    };

    EXPECT_TRUE(detail::MeetsBuchbergerProductCriterion(x_cubed_plus_y, z_fourth));
    const auto s_polynomial = ComputeSPolynomial(x_cubed_plus_y, z_fourth);
    EXPECT_EQ(s_polynomial, expected_s_polynomial);
    EXPECT_EQ(s_polynomial, quotient * y_z_plus_y + expected_fixed_remainder);
    EXPECT_EQ(ReducePolynomial(s_polynomial, ordered_divisors), expected_fixed_remainder);
    EXPECT_TRUE(
        ReducePolynomial(s_polynomial, GradedBasis{z_fourth, y_z_plus_y, x_cubed_plus_y}).IsZero()
    );
}

TEST(BuchbergerAlgorithmTest, LcmRepresentationDoesNotImplyReductionToZero) {
    const PolynomialType f1{
        TermType{Coefficient{1}, Monomial{1, 0, 1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType f2{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType f3{
        TermType{Coefficient{1}, Monomial{1, 0, 1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType expected_s_polynomial{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
    };
    const Basis generators{f1, f2, f3};
    const Monomial leading_lcm =
        ComputeLeastCommonMultiple(f1.GetLeadingMonomial(), f2.GetLeadingMonomial());

    EXPECT_EQ(leading_lcm, (Monomial{1, 1, 1}));
    const auto s_polynomial = ComputeSPolynomial(f1, f2);
    EXPECT_EQ(s_polynomial, expected_s_polynomial);
    EXPECT_EQ(s_polynomial, -f1 + f3);
    EXPECT_TRUE(f1.GetLeadingMonomial().Divides(leading_lcm));
    EXPECT_TRUE(f3.GetLeadingMonomial().Divides(leading_lcm));

    EXPECT_EQ(ReducePolynomial(s_polynomial, generators), s_polynomial);
    const auto completed_basis = ComputeGroebnerBasis(generators);
    EXPECT_GT(completed_basis.size(), generators.size());
    ExpectGroebnerBasis(completed_basis);
}

TEST(BuchbergerAlgorithmTest, TwistedCubicGeneratorsFormAGroebnerBasisWhenZIsHighest) {
    // Degree coordinates are (z, y, x), so lexicographic order is z > y > x.
    const PolynomialType y_minus_x_squared{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
    };
    const PolynomialType z_minus_x_cubed{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 3}},
    };
    const PolynomialType expected_s_polynomial{
        TermType{Coefficient{-1}, Monomial{1, 0, 2}},
        TermType{Coefficient{1}, Monomial{0, 1, 3}},
    };
    const Basis generators{y_minus_x_squared, z_minus_x_cubed};
    const Basis expected_reduced{z_minus_x_cubed, y_minus_x_squared};

    EXPECT_EQ(y_minus_x_squared.GetLeadingMonomial(), (Monomial{0, 1}));
    EXPECT_EQ(z_minus_x_cubed.GetLeadingMonomial(), (Monomial{1}));
    EXPECT_TRUE(detail::MeetsBuchbergerProductCriterion(y_minus_x_squared, z_minus_x_cubed));

    const auto s_polynomial = ComputeSPolynomial(y_minus_x_squared, z_minus_x_cubed);
    EXPECT_EQ(s_polynomial, expected_s_polynomial);
    EXPECT_TRUE(ReducePolynomial(s_polynomial, generators).IsZero());
    EXPECT_EQ(ComputeGroebnerBasis(generators), generators);
    EXPECT_EQ(ComputeReducedGroebnerBasis(generators), expected_reduced);
    ExpectGroebnerBasis(generators);
}

#ifndef NDEBUG
TEST(ComputeSPolynomialDeathTest, ZeroOperandsAreRejected) {
    const PolynomialType nonzero{TermType{Coefficient{1}}};

    EXPECT_DEATH((void)ComputeSPolynomial(PolynomialType{}, nonzero), "nonzero");
    EXPECT_DEATH((void)ComputeSPolynomial(nonzero, PolynomialType{}), "nonzero");
}
#endif

TEST(BuchbergerAlgorithmTest, TextbookSequenceAddsNonzeroSPolynomialRemainders) {
    using GradedPolynomial = Polynomial<Coefficient, GradedLexOrder>;
    using GradedBasis = PolynomialBasis<Coefficient, GradedLexOrder>;

    const GradedPolynomial f1{
        TermType{Coefficient{1}, Monomial{3}},
        TermType{Coefficient{-2}, Monomial{1, 1}},
    };
    const GradedPolynomial f2{
        TermType{Coefficient{1}, Monomial{2, 1}},
        TermType{Coefficient{-2}, Monomial{0, 2}},
        TermType{Coefficient{1}, Monomial{1}},
    };
    const GradedPolynomial x{TermType{Coefficient{1}, Monomial{1}}};
    const GradedPolynomial y{TermType{Coefficient{1}, Monomial{0, 1}}};
    const GradedPolynomial x_squared{TermType{Coefficient{1}, Monomial{2}}};
    const GradedBasis generators{f1, f2};
    const GradedPolynomial f3{TermType{Coefficient{-1}, Monomial{2}}};
    const GradedPolynomial f4{TermType{Coefficient{-2}, Monomial{1, 1}}};
    const GradedPolynomial f5{
        TermType{Coefficient{-2}, Monomial{0, 2}},
        TermType{Coefficient{1}, Monomial{1}},
    };
    const GradedBasis basis_through_f3{f1, f2, f3};
    const GradedBasis basis_through_f4{f1, f2, f3, f4};
    const GradedBasis expected_basis{f1, f2, f3, f4, f5};
    const GradedBasis expected_minimal_basis{
        GradedPolynomial{TermType{Coefficient{1}, Monomial{2}}},
        GradedPolynomial{TermType{Coefficient{1}, Monomial{1, 1}}},
        GradedPolynomial{
            TermType{Coefficient{1}, Monomial{0, 2}},
            TermType{Coefficient{-1, 2}, Monomial{1}},
        },
    };

    EXPECT_EQ(f1.GetLeadingMonomial(), (Monomial{3}));
    EXPECT_EQ(f2.GetLeadingMonomial(), (Monomial{2, 1}));
    EXPECT_FALSE(f1.GetLeadingMonomial().Divides(x_squared.GetLeadingMonomial()));
    EXPECT_FALSE(f2.GetLeadingMonomial().Divides(x_squared.GetLeadingMonomial()));

    EXPECT_EQ(x * f2 - y * f1, x_squared);
    EXPECT_EQ(ComputeSPolynomial(f1, f2), -x_squared);
    EXPECT_EQ(ComputeSPolynomial(f1, f2), f3);
    EXPECT_EQ(ReducePolynomial(ComputeSPolynomial(f1, f2), generators), f3);
    EXPECT_TRUE(ReducePolynomial(ComputeSPolynomial(f1, f2), basis_through_f3).IsZero());

    EXPECT_EQ(ComputeSPolynomial(f1, f3), f4);
    EXPECT_EQ(ReducePolynomial(ComputeSPolynomial(f1, f3), basis_through_f3), f4);
    EXPECT_TRUE(ReducePolynomial(ComputeSPolynomial(f1, f3), basis_through_f4).IsZero());

    EXPECT_EQ(ComputeSPolynomial(f1, f4), y * f4);
    EXPECT_TRUE(ReducePolynomial(ComputeSPolynomial(f1, f4), basis_through_f4).IsZero());

    EXPECT_EQ(ComputeSPolynomial(f2, f3), f5);
    EXPECT_EQ(ReducePolynomial(ComputeSPolynomial(f2, f3), basis_through_f4), f5);
    EXPECT_TRUE(ReducePolynomial(ComputeSPolynomial(f2, f3), expected_basis).IsZero());

    const auto basis = ComputeGroebnerBasis(generators);
    EXPECT_EQ(basis, expected_basis);
    ExpectGroebnerBasis(expected_basis);
    EXPECT_TRUE(ReducePolynomial(x_squared, expected_basis).IsZero());
    EXPECT_TRUE(f3.GetLeadingMonomial().Divides(f1.GetLeadingMonomial()));
    EXPECT_TRUE(f4.GetLeadingMonomial().Divides(f2.GetLeadingMonomial()));
    EXPECT_EQ(detail::MinimizeGroebnerBasis(expected_basis), expected_minimal_basis);
    EXPECT_EQ(ComputeMinimalGroebnerBasis(generators), expected_minimal_basis);
    EXPECT_EQ(ComputeReducedGroebnerBasis(generators), expected_minimal_basis);
}

TEST(BuchbergerAlgorithmTest, LinearGeneratorsWithCoprimeLeadingTermsFormAGroebnerBasis) {
    const PolynomialType x_plus_z{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
    };
    const PolynomialType y_minus_z{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
    };
    const PolynomialType expected_s_polynomial{
        TermType{Coefficient{1}, Monomial{1, 0, 1}},
        TermType{Coefficient{1}, Monomial{0, 1, 1}},
    };
    const Basis generators{x_plus_z, y_minus_z};

    EXPECT_EQ(x_plus_z.GetLeadingMonomial(), (Monomial{1}));
    EXPECT_EQ(y_minus_z.GetLeadingMonomial(), (Monomial{0, 1}));
    EXPECT_TRUE(detail::MeetsBuchbergerProductCriterion(x_plus_z, y_minus_z));

    const auto s_polynomial = ComputeSPolynomial(x_plus_z, y_minus_z);
    EXPECT_EQ(s_polynomial, expected_s_polynomial);
    EXPECT_TRUE(ReducePolynomial(s_polynomial, generators).IsZero());

    EXPECT_EQ(ComputeGroebnerBasis(generators), generators);
    EXPECT_EQ(ComputeReducedGroebnerBasis(generators), generators);
    ExpectGroebnerBasis(generators);
}

TEST(BuchbergerAlgorithmTest, EqualLeadingMonomialsRequireTheDifferenceGenerator) {
    const PolynomialType f1{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
    };
    const PolynomialType f2{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
    };
    const PolynomialType f3{
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
    };
    const Basis original_generators{f1, f2};
    const Basis completed_basis{f1, f2, f3};
    const Basis expected_reduced{
        PolynomialType{
            TermType{Coefficient{1}, Monomial{2}},
            TermType{Coefficient{-1}, Monomial{0, 0, 1}},
        },
        PolynomialType{
            TermType{Coefficient{1}, Monomial{0, 1}},
            TermType{Coefficient{-1}, Monomial{0, 0, 1}},
        },
    };

    EXPECT_EQ(f3, f1 - f2);
    EXPECT_EQ(f1.GetLeadingMonomial(), (Monomial{2}));
    EXPECT_EQ(f2.GetLeadingMonomial(), (Monomial{2}));
    EXPECT_FALSE(f1.GetLeadingMonomial().Divides(f3.GetLeadingMonomial()));
    EXPECT_FALSE(f2.GetLeadingMonomial().Divides(f3.GetLeadingMonomial()));

    EXPECT_EQ(ComputeSPolynomial(f1, f2), f3);
    EXPECT_EQ(ReducePolynomial(f3, original_generators), f3);
    EXPECT_EQ(ComputeGroebnerBasis(original_generators), completed_basis);
    ExpectGroebnerBasis(completed_basis);
    EXPECT_EQ(ComputeReducedGroebnerBasis(original_generators), expected_reduced);
}

TEST(BuchbergerAlgorithmTest, TextbookOverlapReductionAddsTheFourthGenerator) {
    const PolynomialType f1{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 2}},
    };
    const PolynomialType f2{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
    };
    const PolynomialType f3{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 2}},
    };
    const PolynomialType first_s_polynomial{
        TermType{Coefficient{-1}, Monomial{1, 1}},
        TermType{Coefficient{1}, Monomial{1, 0, 1}},
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 2}},
    };
    const PolynomialType after_f2_reductions{
        TermType{Coefficient{2}, Monomial{0, 2}},
        TermType{Coefficient{-2}, Monomial{0, 1, 1}},
        TermType{Coefficient{2}, Monomial{0, 0, 2}},
    };
    const PolynomialType raw_f4{
        TermType{Coefficient{2}, Monomial{0, 0, 4}},
        TermType{Coefficient{2}, Monomial{0, 0, 3}},
        TermType{Coefficient{2}, Monomial{0, 0, 2}},
    };
    const PolynomialType f4{
        TermType{Coefficient{1}, Monomial{0, 0, 4}},
        TermType{Coefficient{1}, Monomial{0, 0, 3}},
        TermType{Coefficient{1}, Monomial{0, 0, 2}},
    };
    const Basis generators{f1, f2, f3};
    const Basis completed_basis{f1, f2, f3, raw_f4};
    const Basis expected_minimal{f2, f3, f4};
    const Basis expected_reduced{
        PolynomialType{
            TermType{Coefficient{1}, Monomial{1}},
            TermType{Coefficient{-1}, Monomial{0, 0, 2}},
            TermType{Coefficient{-1}, Monomial{0, 0, 1}},
        },
        f3,
        f4,
    };

    EXPECT_EQ(ComputeSPolynomial(f1, f2), first_s_polynomial);
    EXPECT_EQ(ReducePolynomial(first_s_polynomial, Basis{f2}), after_f2_reductions);
    EXPECT_EQ(ReducePolynomial(after_f2_reductions, Basis{f3}), raw_f4);
    EXPECT_EQ(ReducePolynomial(first_s_polynomial, Basis{f2, f3}), raw_f4);

    auto normalized_f4 = raw_f4;
    normalized_f4.MakeMonic();
    EXPECT_EQ(normalized_f4, f4);

    EXPECT_TRUE(detail::MeetsBuchbergerProductCriterion(f1, f3));
    EXPECT_TRUE(detail::MeetsBuchbergerProductCriterion(f2, f3));
    EXPECT_TRUE(detail::MeetsBuchbergerProductCriterion(f1, f4));
    EXPECT_TRUE(detail::MeetsBuchbergerProductCriterion(f2, f4));
    EXPECT_TRUE(detail::MeetsBuchbergerProductCriterion(f3, f4));

    EXPECT_EQ(ComputeGroebnerBasis(generators), completed_basis);
    ExpectGroebnerBasis(completed_basis);
    EXPECT_TRUE(f2.GetLeadingMonomial().Divides(f1.GetLeadingMonomial()));
    EXPECT_EQ(detail::MinimizeGroebnerBasis(completed_basis), expected_minimal);
    EXPECT_EQ(ComputeMinimalGroebnerBasis(generators), expected_minimal);
    EXPECT_EQ(ReducePolynomial(f2, Basis{f3, f4}), expected_reduced.front());
    EXPECT_EQ(detail::InterreduceGroebnerBasis(expected_minimal), expected_reduced);
    EXPECT_EQ(ComputeReducedGroebnerBasis(generators), expected_reduced);
}

TEST(BuchbergerAlgorithmTest, EveryImprovementModeProducesAnEquivalentGroebnerBasis) {
    const PolynomialType x_squared_minus_y{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
    };
    const PolynomialType x_y_minus_one{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}},
    };
    const Basis generators{x_squared_minus_y, x_y_minus_one};
    const Basis expected_reduced{
        PolynomialType{
            TermType{Coefficient{1}, Monomial{1}},
            TermType{Coefficient{-1}, Monomial{0, 2}},
        },
        PolynomialType{
            TermType{Coefficient{1}, Monomial{0, 3}},
            TermType{Coefficient{-1}},
        },
    };
    const std::array improvements{
        detail::Improvement::kNone,
        detail::Improvement::kProductCriterion,
        detail::Improvement::kProductAndChainCriteria,
    };

    for (const auto improvement : improvements) {
        SCOPED_TRACE(static_cast<int>(improvement));
        const auto basis = detail::RunBuchbergerAlgorithm(generators, improvement);
        ExpectGroebnerBasis(basis);

        const auto reduced = detail::InterreduceGroebnerBasis(detail::MinimizeGroebnerBasis(basis));
        EXPECT_EQ(reduced, expected_reduced);
    }
}

TEST(BuchbergerAlgorithmTest, ZeroGeneratorsAreRemovedBeforePairsAreBuilt) {
    const PolynomialType generator{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };

    const auto basis =
        detail::RunBuchbergerAlgorithm(Basis{PolynomialType{}, generator, PolynomialType{}});

    ASSERT_EQ(basis.size(), 1U);
    EXPECT_EQ(basis.front(), generator);
}

TEST(MinimizeGroebnerBasisTest, MakesPolynomialsMonicAndRemovesDivisibleLeadingTerms) {
    const Basis basis{
        PolynomialType{TermType{Coefficient{2}, Monomial{2}}},
        PolynomialType{TermType{Coefficient{-3}, Monomial{1}}},
        PolynomialType{TermType{Coefficient{4}, Monomial{0, 1}}},
    };
    const Basis expected{
        PolynomialType{TermType{Coefficient{1}, Monomial{1}}},
        PolynomialType{TermType{Coefficient{1}, Monomial{0, 1}}},
    };

    EXPECT_EQ(detail::MinimizeGroebnerBasis(basis), expected);
}

TEST(MinimizeGroebnerBasisTest, KeepsOnlyTheFirstOfEqualLeadingMonomials) {
    const PolynomialType first{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType second{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{-1}},
    };

    const auto minimal = detail::MinimizeGroebnerBasis(Basis{first, second});

    ASSERT_EQ(minimal.size(), 1U);
    EXPECT_EQ(minimal.front(), first);
}

TEST(MinimizeGroebnerBasisTest, AConstantGeneratorCollapsesTheBasisToOne) {
    const PolynomialType polynomial{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-1}},
    };

    const auto minimal =
        detail::MinimizeGroebnerBasis(Basis{polynomial, PolynomialType{TermType{Coefficient{-7}}}});

    ASSERT_EQ(minimal.size(), 1U);
    EXPECT_TRUE(minimal.front().IsOne());
}

TEST(GroebnerBasisDetailTest, SortOrdersByDescendingLeadingMonomial) {
    Basis basis{
        PolynomialType{TermType{Coefficient{1}}},
        PolynomialType{TermType{Coefficient{1}, Monomial{0, 1}}},
        PolynomialType{TermType{Coefficient{1}, Monomial{2}}},
    };
    const Basis expected{
        PolynomialType{TermType{Coefficient{1}, Monomial{2}}},
        PolynomialType{TermType{Coefficient{1}, Monomial{0, 1}}},
        PolynomialType{TermType{Coefficient{1}}},
    };

    detail::SortBasisByLeadingMonomial(basis);

    EXPECT_EQ(basis, expected);
}

TEST(GroebnerBasisDetailTest, InterreduceRemovesTermsReducibleByOtherBasisElements) {
    const PolynomialType x_plus_y{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };
    const PolynomialType y{TermType{Coefficient{1}, Monomial{0, 1}}};
    const Basis expected{
        PolynomialType{TermType{Coefficient{1}, Monomial{1}}},
        y,
    };

    EXPECT_EQ(detail::InterreduceGroebnerBasis(Basis{x_plus_y, y}), expected);
}

TEST(ComputeGroebnerBasisTest, ResultPassesBuchbergersCriterion) {
    const PolynomialType first{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
    };
    const PolynomialType second{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}},
    };

    const auto basis = ComputeGroebnerBasis(Basis{first, second});

    EXPECT_GE(basis.size(), 2U);
    ExpectGroebnerBasis(basis);
    EXPECT_TRUE(ReducePolynomial(first, basis).IsZero());
    EXPECT_TRUE(ReducePolynomial(second, basis).IsZero());
}

TEST(ComputeGroebnerBasisTest, SymmetricSystemMatchesTheTextbookLexBasis) {
    const PolynomialType f1{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType f2{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType f3{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 2}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType g1 = f3;
    const PolynomialType g2{
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
    };
    const PolynomialType g3{
        TermType{Coefficient{2}, Monomial{0, 1, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 4}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
    };
    const PolynomialType g4{
        TermType{Coefficient{1}, Monomial{0, 0, 6}},
        TermType{Coefficient{-4}, Monomial{0, 0, 4}},
        TermType{Coefficient{4}, Monomial{0, 0, 3}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
    };
    const PolynomialType reduced_g3{
        TermType{Coefficient{1}, Monomial{0, 1, 2}},
        TermType{Coefficient{1, 2}, Monomial{0, 0, 4}},
        TermType{Coefficient{-1, 2}, Monomial{0, 0, 2}},
    };
    const Basis generators{f1, f2, f3};
    const Basis textbook_basis{g1, g2, g3, g4};
    const Basis expected_reduced{g1, g2, reduced_g3, g4};

    ExpectGroebnerBasis(textbook_basis);
    for (const auto& generator : generators) {
        EXPECT_TRUE(ReducePolynomial(generator, textbook_basis).IsZero());
    }

    EXPECT_EQ(ComputeReducedGroebnerBasis(generators), expected_reduced);
    EXPECT_EQ(ComputeReducedGroebnerBasis(textbook_basis), expected_reduced);

    const PolynomialType z_squared{TermType{Coefficient{1}, Monomial{0, 0, 2}}};
    const PolynomialType z_minus_one{
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType final_factor{
        TermType{Coefficient{1}, Monomial{0, 0, 2}},
        TermType{Coefficient{2}, Monomial{0, 0, 1}},
        TermType{Coefficient{-1}},
    };
    EXPECT_EQ(g4, z_squared * z_minus_one * z_minus_one * final_factor);
}

TEST(ComputeGroebnerBasisTest, TangentSurfaceEliminationMatchesTheTextbookLexBasis) {
    // Degree coordinates are (t, u, x, y, z).
    const PolynomialType first_generator{
        TermType{Coefficient{-1}, Monomial{1}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
    };
    const PolynomialType second_generator{
        TermType{Coefficient{-1}, Monomial{2}},
        TermType{Coefficient{-2}, Monomial{1, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 0, 1}},
    };
    const PolynomialType third_generator{
        TermType{Coefficient{-1}, Monomial{3}},
        TermType{Coefficient{-3}, Monomial{2, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 0, 0, 1}},
    };
    const PolynomialType g1{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
    };
    const PolynomialType g2{
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 0, 1}},
    };
    const PolynomialType g3{
        TermType{Coefficient{1}, Monomial{0, 1, 2}},
        TermType{Coefficient{-1}, Monomial{0, 1, 0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 3}},
        TermType{Coefficient{3, 2}, Monomial{0, 0, 1, 1}},
        TermType{Coefficient{-1, 2}, Monomial{0, 0, 0, 0, 1}},
    };
    const PolynomialType g4{
        TermType{Coefficient{1}, Monomial{0, 1, 1, 1}},
        TermType{Coefficient{-1}, Monomial{0, 1, 0, 0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1, 0, 1}},
        TermType{Coefficient{2}, Monomial{0, 0, 0, 2}},
    };
    const PolynomialType g5{
        TermType{Coefficient{1}, Monomial{0, 1, 1, 0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 1, 0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 2, 0, 1}},
        TermType{Coefficient{-1, 2}, Monomial{0, 0, 1, 2}},
        TermType{Coefficient{-1, 2}, Monomial{0, 0, 0, 1, 1}},
    };
    const PolynomialType g6{
        TermType{Coefficient{1}, Monomial{0, 1, 0, 3}},
        TermType{Coefficient{-1}, Monomial{0, 1, 0, 0, 2}},
        TermType{Coefficient{-2}, Monomial{0, 0, 2, 1, 1}},
        TermType{Coefficient{1, 2}, Monomial{0, 0, 1, 3}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1, 0, 2}},
        TermType{Coefficient{5, 2}, Monomial{0, 0, 0, 2, 1}},
    };
    const PolynomialType g7{
        TermType{Coefficient{1}, Monomial{0, 0, 3, 0, 1}},
        TermType{Coefficient{-3, 4}, Monomial{0, 0, 2, 2}},
        TermType{Coefficient{-3, 2}, Monomial{0, 0, 1, 1, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 0, 3}},
        TermType{Coefficient{1, 4}, Monomial{0, 0, 0, 0, 2}},
    };
    const Basis generators{
        first_generator,
        second_generator,
        third_generator,
    };
    const Basis expected_basis{g1, g2, g3, g4, g5, g6, g7};

    ExpectGroebnerBasis(expected_basis);
    for (const auto& generator : generators) {
        EXPECT_TRUE(ReducePolynomial(generator, expected_basis).IsZero());
    }
    EXPECT_EQ(ComputeReducedGroebnerBasis(expected_basis), expected_basis);
    EXPECT_EQ(ComputeReducedGroebnerBasis(generators), expected_basis);

    EXPECT_EQ(g7.GetLeadingMonomial().GetDegreeOf(0), 0U);
    EXPECT_EQ(g7.GetLeadingMonomial().GetDegreeOf(1), 0U);
    for (const auto& term : g7.GetTerms()) {
        EXPECT_EQ(term.GetMonomial().GetDegreeOf(0), 0U);
        EXPECT_EQ(term.GetMonomial().GetDegreeOf(1), 0U);
    }
}

TEST(ComputeGroebnerBasisTest, FirstEliminationIdealMatchesTheTextbookLexBasis) {
    const PolynomialType first_generator{
        TermType{Coefficient{1}, Monomial{2, 1}},
        TermType{Coefficient{1}, Monomial{1, 0, 1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType second_generator{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}, Monomial{1, 0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType g1{
        TermType{Coefficient{1}, Monomial{2, 0, 2}},
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const PolynomialType g2 = second_generator;
    const PolynomialType g3{
        TermType{Coefficient{1}, Monomial{1, 0, 3}}, TermType{Coefficient{-1}, Monomial{1, 0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 1}},   TermType{Coefficient{1}, Monomial{0, 0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}}, TermType{Coefficient{-1}},
    };
    const PolynomialType g4{
        TermType{Coefficient{1}, Monomial{0, 2}},     TermType{Coefficient{-2}, Monomial{0, 1, 2}},
        TermType{Coefficient{-1}, Monomial{0, 1, 1}}, TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{2}, Monomial{0, 0, 4}},  TermType{Coefficient{-1}, Monomial{0, 0, 3}},
    };
    const Basis generators{first_generator, second_generator};
    const Basis expected_basis{g1, g2, g3, g4};

    ExpectGroebnerBasis(expected_basis);
    for (const auto& generator : generators) {
        EXPECT_TRUE(ReducePolynomial(generator, expected_basis).IsZero());
    }
    EXPECT_EQ(ComputeReducedGroebnerBasis(generators), expected_basis);

    const PolynomialType x{TermType{Coefficient{1}, Monomial{1}}};
    const PolynomialType y_minus_z_squared{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
    };
    const PolynomialType z_cubed_minus_z_squared{
        TermType{Coefficient{1}, Monomial{0, 0, 3}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
    };
    const PolynomialType g3_tail{
        TermType{Coefficient{-1}, Monomial{0, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType z_minus_one{
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{-1}},
    };
    EXPECT_EQ(g2, y_minus_z_squared * x + z_minus_one);
    EXPECT_EQ(g3, z_cubed_minus_z_squared * x + g3_tail);

    EXPECT_GT(g1.GetLeadingMonomial().GetDegreeOf(0), 0U);
    EXPECT_GT(g2.GetLeadingMonomial().GetDegreeOf(0), 0U);
    EXPECT_GT(g3.GetLeadingMonomial().GetDegreeOf(0), 0U);
    for (const auto& term : g4.GetTerms()) {
        EXPECT_EQ(term.GetMonomial().GetDegreeOf(0), 0U);
    }
    EXPECT_EQ(ComputeReducedGroebnerBasis(Basis{g4}), Basis{g4});
}

TEST(ComputeGroebnerBasisTest, ThreeEquationSystemHasTheExpectedSolutions) {
    // Degree coordinates are (a, b, c).
    const PolynomialType f1{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
    };
    const PolynomialType f2{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-1}, Monomial{1}},
        TermType{Coefficient{-1}, Monomial{0, 1, 1}},
    };
    const PolynomialType f3{
        TermType{Coefficient{1}, Monomial{1, 0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
    };
    const PolynomialType f4{
        TermType{Coefficient{-1}, Monomial{0, 2}},    TermType{Coefficient{-1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 1, 1}}, TermType{Coefficient{1}, Monomial{0, 0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
    };
    const PolynomialType f5{TermType{Coefficient{2}, Monomial{0, 1, 1}}};
    const PolynomialType f6{
        TermType{Coefficient{-2}, Monomial{0, 0, 3}},
        TermType{Coefficient{-2}, Monomial{0, 0, 2}},
    };
    const Basis generators{f1, f2, f3};
    const Basis completed_basis{f1, f2, f3, f4, f5, f6};
    const Basis expected_reduced{
        PolynomialType{
            TermType{Coefficient{1}, Monomial{2}},
            TermType{Coefficient{-1}, Monomial{1}},
        },
        f1,
        PolynomialType{
            TermType{Coefficient{1}, Monomial{1, 0, 1}},
            TermType{Coefficient{-1}, Monomial{0, 0, 2}},
            TermType{Coefficient{-1}, Monomial{0, 0, 1}},
        },
        PolynomialType{
            TermType{Coefficient{1}, Monomial{0, 2}},
            TermType{Coefficient{1}, Monomial{0, 1}},
            TermType{Coefficient{-1}, Monomial{0, 0, 2}},
            TermType{Coefficient{-1}, Monomial{0, 0, 1}},
        },
        PolynomialType{TermType{Coefficient{1}, Monomial{0, 1, 1}}},
        PolynomialType{
            TermType{Coefficient{1}, Monomial{0, 0, 3}},
            TermType{Coefficient{1}, Monomial{0, 0, 2}},
        },
    };

    EXPECT_EQ(ReducePolynomial(ComputeSPolynomial(f1, f2), Basis{f1, f2, f3}), f4);
    EXPECT_EQ(ReducePolynomial(ComputeSPolynomial(f2, f3), Basis{f1, f2, f3, f4}), f5);
    EXPECT_EQ(ReducePolynomial(ComputeSPolynomial(f1, f3), Basis{f1, f2, f3, f4, f5}), f6);

    ExpectGroebnerBasis(completed_basis);
    EXPECT_EQ(ComputeReducedGroebnerBasis(generators), expected_reduced);

    const PolynomialType c_squared{TermType{Coefficient{1}, Monomial{0, 0, 2}}};
    const PolynomialType c_plus_one{
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{1}},
    };
    EXPECT_EQ(expected_reduced.back(), c_squared * c_plus_one);

    const auto is_solution = [](int a_value, int b_value, int c_value) {
        const Coefficient a{a_value};
        const Coefficient b{b_value};
        const Coefficient c{c_value};
        return a * b == c * c + c && a * a == a + b * c && a * c == b * b + b;
    };

    std::vector<std::array<int, 3>> solutions;
    for (const int c_value : {0, -1}) {
        for (const int a_value : {0, 1}) {
            for (const int b_value : {0, -1}) {
                if (is_solution(a_value, b_value, c_value)) {
                    solutions.push_back({a_value, b_value, c_value});
                }
            }
        }
    }
    std::sort(solutions.begin(), solutions.end());

    const std::vector<std::array<int, 3>> expected_solutions{
        {0, -1, 0},
        {0, 0, -1},
        {0, 0, 0},
        {1, 0, 0},
    };
    EXPECT_EQ(solutions, expected_solutions);
    EXPECT_FALSE(is_solution(1, 0, -1));
    EXPECT_FALSE(is_solution(1, -1, 0));
}

TEST(ComputeGroebnerBasisTest, PositiveLinearATermLeavesAFreeVariable) {
    // Degree coordinates are (a, b, c).
    const PolynomialType f1{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
    };
    const PolynomialType f2{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{-1}, Monomial{0, 1, 1}},
    };
    const PolynomialType f3{
        TermType{Coefficient{1}, Monomial{1, 0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
    };
    const PolynomialType f4{
        TermType{Coefficient{-1}, Monomial{0, 2}},    TermType{Coefficient{-1}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 1, 1}}, TermType{Coefficient{-1}, Monomial{0, 0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
    };
    const Basis generators{f1, f2, f3};
    const Basis completed_basis{f1, f2, f3, f4};
    const PolynomialType monic_f4 = -f4;
    const Basis expected_minimal{f1, f2, f3, monic_f4};
    const Basis expected_reduced{
        f2,
        f1,
        PolynomialType{
            TermType{Coefficient{1}, Monomial{1, 0, 1}},
            TermType{Coefficient{1}, Monomial{0, 1, 1}},
            TermType{Coefficient{1}, Monomial{0, 0, 2}},
            TermType{Coefficient{1}, Monomial{0, 0, 1}},
        },
        PolynomialType{
            TermType{Coefficient{1}, Monomial{0, 2}},
            TermType{Coefficient{1}, Monomial{0, 1}},
            TermType{Coefficient{1}, Monomial{0, 1, 1}},
            TermType{Coefficient{1}, Monomial{0, 0, 2}},
            TermType{Coefficient{1}, Monomial{0, 0, 1}},
        },
    };

    EXPECT_EQ(ReducePolynomial(ComputeSPolynomial(f1, f2), generators), f4);
    EXPECT_EQ(ComputeGroebnerBasis(generators), completed_basis);
    ExpectGroebnerBasis(completed_basis);
    EXPECT_EQ(detail::MinimizeGroebnerBasis(completed_basis), expected_minimal);
    EXPECT_EQ(ComputeMinimalGroebnerBasis(generators), expected_minimal);
    for (std::size_t lhs = 0; lhs < expected_minimal.size(); ++lhs) {
        EXPECT_TRUE(expected_minimal[lhs].IsMonic());
        for (std::size_t rhs = 0; rhs < expected_minimal.size(); ++rhs) {
            if (lhs != rhs) {
                EXPECT_FALSE(
                    expected_minimal[lhs].GetLeadingMonomial().Divides(
                        expected_minimal[rhs].GetLeadingMonomial()
                    )
                );
            }
        }
    }
    EXPECT_NE(ReducePolynomial(f3, Basis{monic_f4}), f3);
    EXPECT_EQ(ReducePolynomial(f3, Basis{monic_f4}), expected_reduced[2]);
    EXPECT_EQ(detail::InterreduceGroebnerBasis(expected_minimal), expected_reduced);
    EXPECT_EQ(ComputeReducedGroebnerBasis(generators), expected_reduced);

    for (const auto& polynomial : expected_reduced) {
        const auto terms = polynomial.GetTerms();
        const bool contains_only_c = std::all_of(terms.begin(), terms.end(), [](const auto& term) {
            return term.GetMonomial().GetDegreeOf(0) == 0 && term.GetMonomial().GetDegreeOf(1) == 0;
        });
        EXPECT_FALSE(contains_only_c);
    }
    EXPECT_TRUE(
        std::any_of(expected_reduced.begin(), expected_reduced.end(), [](const auto& polynomial) {
            const auto terms = polynomial.GetTerms();
            return std::all_of(terms.begin(), terms.end(), [](const auto& term) {
                return term.GetMonomial().GetDegreeOf(0) == 0;
            });
        })
    );

    const auto is_solution = [](const Coefficient& a, const Coefficient& b, const Coefficient& c) {
        return a * b == c * c + c && a * a + a == b * c && a * c == b * b + b;
    };
    EXPECT_TRUE(is_solution(Coefficient{0}, Coefficient{0}, Coefficient{0}));
    EXPECT_TRUE(is_solution(Coefficient{-1}, Coefficient{0}, Coefficient{0}));
    EXPECT_TRUE(is_solution(Coefficient{0}, Coefficient{-1}, Coefficient{0}));
    EXPECT_TRUE(is_solution(Coefficient{0}, Coefficient{0}, Coefficient{-1}));
    EXPECT_TRUE(is_solution(Coefficient{-2, 3}, Coefficient{-2, 3}, Coefficient{1, 3}));
}

TEST(ComputeGroebnerBasisTest, VariablePriorityRevealsDifferentDependentPairs) {
    const auto contains_polynomial_without_first_variable = [](const auto& basis) {
        return std::any_of(basis.begin(), basis.end(), [](const auto& polynomial) {
            const auto terms = polynomial.GetTerms();
            return std::all_of(terms.begin(), terms.end(), [](const auto& term) {
                return term.GetMonomial().GetDegreeOf(0) == 0;
            });
        });
    };

    // Degree coordinates are (b, a, c), so lexicographic order is b > a > c.
    const PolynomialType ba_f1{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
    };
    const PolynomialType ba_f2{
        TermType{Coefficient{1}, Monomial{1, 0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
    };
    const PolynomialType ba_f3{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{-1}, Monomial{0, 1, 1}},
    };
    const PolynomialType a_c_relation{
        TermType{Coefficient{1}, Monomial{0, 3}},
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 0, 3}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
    };
    const Basis ba_generators{ba_f1, ba_f2, ba_f3};

    EXPECT_EQ(ComputeSPolynomial(ba_f1, ba_f2), a_c_relation);
    EXPECT_EQ(ReducePolynomial(a_c_relation, ba_generators), a_c_relation);
    const auto ba_basis = ComputeGroebnerBasis(ba_generators);
    ExpectGroebnerBasis(ba_basis);
    EXPECT_TRUE(ReducePolynomial(a_c_relation, ba_basis).IsZero());
    EXPECT_TRUE(
        contains_polynomial_without_first_variable(ComputeReducedGroebnerBasis(ba_generators))
    );

    // Degree coordinates are (c, a, b), so lexicographic order is c > a > b.
    const PolynomialType ca_f1{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{-1}, Monomial{0, 1, 1}},
    };
    const PolynomialType ca_f2{
        TermType{Coefficient{1}, Monomial{1, 0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
    };
    const PolynomialType ca_f3{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
    };
    const PolynomialType a_b_relation{
        TermType{Coefficient{-1}, Monomial{0, 3}},
        TermType{Coefficient{-1}, Monomial{0, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 3}},
        TermType{Coefficient{1}, Monomial{0, 0, 2}},
    };
    const Basis ca_generators{ca_f1, ca_f2, ca_f3};

    EXPECT_EQ(ComputeSPolynomial(ca_f2, ca_f3), a_b_relation);
    EXPECT_EQ(ReducePolynomial(a_b_relation, ca_generators), a_b_relation);
    const auto ca_basis = ComputeGroebnerBasis(ca_generators);
    ExpectGroebnerBasis(ca_basis);
    EXPECT_TRUE(ReducePolynomial(a_b_relation, ca_basis).IsZero());
    EXPECT_TRUE(
        contains_polynomial_without_first_variable(ComputeReducedGroebnerBasis(ca_generators))
    );
}

TEST(ComputeMinimalGroebnerBasisTest, RemovesRedundantInputGenerators) {
    const PolynomialType x_squared{TermType{Coefficient{2}, Monomial{2}}};
    const PolynomialType x{TermType{Coefficient{-3}, Monomial{1}}};
    const PolynomialType y{TermType{Coefficient{4}, Monomial{0, 1}}};
    const Basis expected{
        PolynomialType{TermType{Coefficient{1}, Monomial{1}}},
        PolynomialType{TermType{Coefficient{1}, Monomial{0, 1}}},
    };

    EXPECT_EQ(ComputeMinimalGroebnerBasis(Basis{x_squared, x, y}), expected);
}

TEST(ComputeMinimalGroebnerBasisTest, MinimalBasisNeedNotBeUnique) {
    using GradedPolynomial = Polynomial<Coefficient, GradedLexOrder>;
    using GradedBasis = PolynomialBasis<Coefficient, GradedLexOrder>;

    const GradedPolynomial x_squared{TermType{Coefficient{1}, Monomial{2}}};
    const GradedPolynomial x_y{TermType{Coefficient{1}, Monomial{1, 1}}};
    const GradedPolynomial y_squared_minus_half_x{
        TermType{Coefficient{1}, Monomial{0, 2}},
        TermType{Coefficient{-1, 2}, Monomial{1}},
    };
    const GradedPolynomial x_squared_plus_three_x_y{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{3}, Monomial{1, 1}},
    };
    const GradedBasis canonical_basis{
        x_squared,
        x_y,
        y_squared_minus_half_x,
    };
    const GradedBasis alternate_basis{
        x_squared_plus_three_x_y,
        x_y,
        y_squared_minus_half_x,
    };

    EXPECT_NE(alternate_basis, canonical_basis);
    ExpectGroebnerBasis(alternate_basis);
    EXPECT_EQ(detail::MinimizeGroebnerBasis(alternate_basis), alternate_basis);

    for (std::size_t lhs = 0; lhs < alternate_basis.size(); ++lhs) {
        EXPECT_TRUE(alternate_basis[lhs].IsMonic());
        for (std::size_t rhs = 0; rhs < alternate_basis.size(); ++rhs) {
            if (lhs != rhs) {
                EXPECT_FALSE(
                    alternate_basis[lhs].GetLeadingMonomial().Divides(
                        alternate_basis[rhs].GetLeadingMonomial()
                    )
                );
            }
        }
    }

    for (const auto& polynomial : canonical_basis) {
        EXPECT_TRUE(ReducePolynomial(polynomial, alternate_basis).IsZero());
    }
    for (const auto& polynomial : alternate_basis) {
        EXPECT_TRUE(ReducePolynomial(polynomial, canonical_basis).IsZero());
    }

    EXPECT_EQ(ComputeReducedGroebnerBasis(alternate_basis), canonical_basis);
}

TEST(ComputeMinimalGroebnerBasisTest, LinearSystemMatchesRowEchelonGroebnerBases) {
    const PolynomialType first_generator{
        TermType{Coefficient{3}, Monomial{1}},
        TermType{Coefficient{-6}, Monomial{0, 1}},
        TermType{Coefficient{-2}, Monomial{0, 0, 1}},
    };
    const PolynomialType second_generator{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{-4}, Monomial{0, 1}},
        TermType{Coefficient{4}, Monomial{0, 0, 0, 1}},
    };
    const PolynomialType third_generator{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{-2}, Monomial{0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 0, 0, 1}},
    };
    const PolynomialType z_plus_three_w{
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{3}, Monomial{0, 0, 0, 1}},
    };
    const PolynomialType reduced_first{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{-2}, Monomial{0, 1}},
        TermType{Coefficient{2}, Monomial{0, 0, 0, 1}},
    };
    const Basis generators{
        first_generator,
        second_generator,
        third_generator,
    };
    const Basis expected_minimal{
        third_generator,
        z_plus_three_w,
    };
    const Basis expected_reduced{
        reduced_first,
        z_plus_three_w,
    };

    EXPECT_EQ(
        z_plus_three_w, (Coefficient{3, 4} * second_generator - Coefficient{1, 2} * first_generator)
    );
    EXPECT_EQ(first_generator, Coefficient{3} * third_generator + z_plus_three_w);
    EXPECT_EQ(second_generator, Coefficient{2} * third_generator + Coefficient{2} * z_plus_three_w);

    ExpectGroebnerBasis(expected_minimal);
    EXPECT_TRUE(expected_minimal[0].IsMonic());
    EXPECT_TRUE(expected_minimal[1].IsMonic());
    EXPECT_FALSE(
        expected_minimal[0].GetLeadingMonomial().Divides(expected_minimal[1].GetLeadingMonomial())
    );
    EXPECT_FALSE(
        expected_minimal[1].GetLeadingMonomial().Divides(expected_minimal[0].GetLeadingMonomial())
    );
    for (const auto& generator : generators) {
        EXPECT_TRUE(ReducePolynomial(generator, expected_minimal).IsZero());
    }

    EXPECT_EQ(
        ComputeMinimalGroebnerBasis(Basis{third_generator, first_generator, second_generator}),
        expected_minimal
    );
    EXPECT_EQ(detail::InterreduceGroebnerBasis(expected_minimal), expected_reduced);
    EXPECT_EQ(ComputeReducedGroebnerBasis(generators), expected_reduced);
}

template <typename Order>
void ExpectCoordinateIdealHasCanonicalReducedBasis() {
    using OrderedPolynomial = Polynomial<Coefficient, Order>;
    using OrderedBasis = PolynomialBasis<Coefficient, Order>;

    const OrderedPolynomial first{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{2}, Monomial{0, 1}},
    };
    const OrderedPolynomial second{
        TermType{Coefficient{-3}, Monomial{1}},
        TermType{Coefficient{3}, Monomial{0, 1}},
    };
    const OrderedBasis expected{
        OrderedPolynomial{TermType{Coefficient{1}, Monomial{1}}},
        OrderedPolynomial{TermType{Coefficient{1}, Monomial{0, 1}}},
    };

    const auto reduced = ComputeReducedGroebnerBasis(OrderedBasis{second, first});

    EXPECT_EQ(reduced, expected);
    for (const auto& polynomial : reduced) {
        EXPECT_TRUE(polynomial.IsMonic());
    }
    ExpectGroebnerBasis(reduced);
}

TEST(ComputeReducedGroebnerBasisTest, IsCanonicalAcrossSupportedMonomialOrders) {
    ExpectCoordinateIdealHasCanonicalReducedBasis<LexicographicOrder>();
    ExpectCoordinateIdealHasCanonicalReducedBasis<GradedLexOrder>();
    ExpectCoordinateIdealHasCanonicalReducedBasis<GradedReverseLexOrder>();
}

TEST(ComputeReducedGroebnerBasisTest, SupportsPrimeFieldCoefficients) {
    using F5 = fields::PrimeField<5>;
    using F5Term = Term<F5>;
    using F5Polynomial = Polynomial<F5, LexicographicOrder>;
    using F5Basis = PolynomialBasis<F5, LexicographicOrder>;

    const F5Polynomial first{
        F5Term{F5{1}, Monomial{1}},
        F5Term{F5{1}, Monomial{0, 1}},
    };
    const F5Polynomial second{
        F5Term{F5{1}, Monomial{1}},
        F5Term{F5{-1}, Monomial{0, 1}},
    };
    const F5Basis expected{
        F5Polynomial{F5Term{F5{1}, Monomial{1}}},
        F5Polynomial{F5Term{F5{1}, Monomial{0, 1}}},
    };

    const auto reduced = ComputeReducedGroebnerBasis(F5Basis{first, second});

    EXPECT_EQ(reduced, expected);
    ExpectGroebnerBasis(reduced);
}

TEST(ComputeReducedGroebnerBasisTest, MoraFamilyContainsTheHighDegreePolynomial) {
    using GrevlexPolynomial = Polynomial<Coefficient, GradedReverseLexOrder>;
    using GrevlexBasis = PolynomialBasis<Coefficient, GradedReverseLexOrder>;

    for (Monomial::Degree n = 1; n <= 3; ++n) {
        SCOPED_TRACE("n = " + std::to_string(n));

        const GrevlexPolynomial first{
            TermType{Coefficient{1}, Monomial{n + 1}},
            TermType{
                Coefficient{-1},
                Monomial{0, 1, n - 1, 1},
            },
        };
        const GrevlexPolynomial second{
            TermType{Coefficient{1}, Monomial{1, n - 1}},
            TermType{Coefficient{-1}, Monomial{0, 0, n}},
        };
        const GrevlexPolynomial third{
            TermType{Coefficient{1}, Monomial{n, 0, 1}},
            TermType{Coefficient{-1}, Monomial{0, n, 0, 1}},
        };
        const GrevlexPolynomial expected{
            TermType{Coefficient{1}, Monomial{0, 0, n * n + 1}},
            TermType{Coefficient{-1}, Monomial{0, n * n, 0, 1}},
        };

        const auto reduced = ComputeReducedGroebnerBasis(GrevlexBasis{first, second, third});

        EXPECT_EQ(expected.GetLeadingMonomial().GetTotalDegree(), n * n + 1);
        EXPECT_NE(std::find(reduced.begin(), reduced.end(), expected), reduced.end());
        EXPECT_TRUE(ReducePolynomial(expected, reduced).IsZero());
        ExpectGroebnerBasis(reduced);
    }
}

TEST(ComputeReducedGroebnerBasisTest, HandlesZeroAndUnitIdeals) {
    EXPECT_TRUE(ComputeReducedGroebnerBasis(Basis{}).empty());
    EXPECT_TRUE(ComputeReducedGroebnerBasis(Basis{PolynomialType{}, PolynomialType{}}).empty());

    const PolynomialType polynomial{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const auto unit_basis =
        ComputeReducedGroebnerBasis(Basis{polynomial, PolynomialType{TermType{Coefficient{5}}}});

    ASSERT_EQ(unit_basis.size(), 1U);
    EXPECT_TRUE(unit_basis.front().IsOne());
}

}  // namespace
}  // namespace groebner_basis::algorithms
