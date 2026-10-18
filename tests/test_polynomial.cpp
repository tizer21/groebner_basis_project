#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include "groebner_basis/fields/rational.h"
#include "groebner_basis/monomial_orders.h"
#include "groebner_basis/polynomial.h"

namespace groebner_basis {
namespace {

using Coefficient = fields::Rational;
using TermType = Term<Coefficient>;
using PolynomialType = Polynomial<Coefficient, LexicographicOrder>;

TEST(PolynomialTest, DefaultConstructionCreatesTheZeroPolynomial) {
    const PolynomialType polynomial;

    EXPECT_TRUE(polynomial.IsZero());
    EXPECT_FALSE(polynomial.IsOne());
    EXPECT_FALSE(polynomial.IsMonic());
    EXPECT_TRUE(polynomial.GetTerms().empty());

    std::ostringstream output;
    output << polynomial;
    EXPECT_EQ(output.str(), "0");
}

#ifndef NDEBUG
TEST(PolynomialDeathTest, ZeroPolynomialHasNoLeadingData) {
    EXPECT_DEATH(
        {
            const PolynomialType zero;
            (void)zero.GetLeadingCoefficient();
        },
        "no leading coefficient"
    );
    EXPECT_DEATH(
        {
            const PolynomialType zero;
            (void)zero.GetLeadingMonomial();
        },
        "no leading monomial"
    );
    EXPECT_DEATH(
        {
            const PolynomialType zero;
            (void)zero.GetLeadingTerm();
        },
        "no leading term"
    );
}
#endif

TEST(PolynomialTest, ConstructionCombinesLikeTermsAndRemovesZeroCoefficients) {
    const PolynomialType polynomial{
        TermType{Coefficient{2}, Monomial{2}},    TermType{Coefficient{-2}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{1, 1}}, TermType{Coefficient{2}, Monomial{1, 1}},
        TermType{Coefficient{0}, Monomial{9, 9}}, TermType{Coefficient{4}},
    };
    const std::vector<TermType> expected{
        TermType{Coefficient{3}, Monomial{1, 1}},
        TermType{Coefficient{4}},
    };

    EXPECT_EQ(polynomial.GetTerms(), expected);
    EXPECT_EQ(polynomial.GetLeadingTerm(), expected.front());
}

TEST(PolynomialTest, VectorConstructionDoesNotModifyItsInput) {
    const std::vector<TermType> terms{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{-1}},
    };
    const auto original = terms;
    const PolynomialType polynomial{terms};

    EXPECT_EQ(terms, original);
    EXPECT_EQ(polynomial.GetTerms(), original);
}

TEST(PolynomialTest, LeadingTermDependsOnTheSelectedMonomialOrder) {
    const std::vector<TermType> terms{
        TermType{Coefficient{3}, Monomial{3, 2}},
        TermType{Coefficient{-5}, Monomial{2, 1, 3}},
    };
    const Polynomial<Coefficient, LexicographicOrder> lexicographic{terms};
    const Polynomial<Coefficient, GradedLexOrder> graded_lexicographic{terms};
    const Polynomial<Coefficient, GradedReverseLexOrder> graded_reverse{terms};

    EXPECT_EQ(lexicographic.GetLeadingTerm(), terms[0]);
    EXPECT_EQ(graded_lexicographic.GetLeadingTerm(), terms[1]);
    EXPECT_EQ(graded_reverse.GetLeadingTerm(), terms[1]);
    EXPECT_EQ(lexicographic.GetLeadingCoefficient(), Coefficient{3});
    EXPECT_EQ(lexicographic.GetLeadingMonomial(), (Monomial{3, 2}));
}

TEST(PolynomialTest, LexLeadingTermPrioritizesTheFirstDifferingVariableOverDegree) {
    const TermType three_x1{Coefficient{3}, Monomial{1}};
    const TermType x2_squared_x3_fourth{
        Coefficient{1},
        Monomial{0, 2, 4},
    };
    const TermType x4_seventh{
        Coefficient{1},
        Monomial{0, 0, 0, 7},
    };
    const std::vector<TermType> first_terms{
        three_x1,
        x2_squared_x3_fourth,
        x4_seventh,
    };
    const PolynomialType first_lexicographic{first_terms};
    const Polynomial<Coefficient, GradedLexOrder> first_graded{first_terms};

    EXPECT_EQ(first_lexicographic.GetLeadingTerm(), three_x1);
    EXPECT_EQ(first_lexicographic.GetLeadingCoefficient(), Coefficient{3});
    EXPECT_EQ(first_lexicographic.GetLeadingMonomial(), (Monomial{1}));
    EXPECT_EQ(first_graded.GetLeadingTerm(), x4_seventh);

    const TermType x_squared_y{Coefficient{1}, Monomial{2, 1}};
    const TermType x_y_squared{Coefficient{1}, Monomial{1, 2}};
    const TermType z_fifth{Coefficient{1}, Monomial{0, 0, 5}};
    const PolynomialType second_lexicographic{
        x_squared_y,
        x_y_squared,
        z_fifth,
    };

    EXPECT_EQ(second_lexicographic.GetLeadingTerm(), x_squared_y);
    EXPECT_EQ(second_lexicographic.GetLeadingMonomial(), (Monomial{2, 1}));
}

TEST(PolynomialTest, TextbookExampleRepresentationFollowsTheSelectedMonomialOrder) {
    const TermType x_y_squared_z{Coefficient{4}, Monomial{1, 2, 1}};
    const TermType z_squared{Coefficient{4}, Monomial{0, 0, 2}};
    const TermType negative_x_cubed{Coefficient{-5}, Monomial{3}};
    const TermType x_squared_z_squared{Coefficient{7}, Monomial{2, 0, 2}};
    const std::vector<TermType> terms{
        x_y_squared_z,
        z_squared,
        negative_x_cubed,
        x_squared_z_squared,
    };

    const Polynomial<Coefficient, LexicographicOrder> lexicographic{terms};
    const Polynomial<Coefficient, GradedLexOrder> graded_lexicographic{terms};
    const Polynomial<Coefficient, GradedReverseLexOrder> graded_reverse{terms};

    EXPECT_EQ(
        lexicographic.GetTerms(), (std::vector<TermType>{
                                      negative_x_cubed,
                                      x_squared_z_squared,
                                      x_y_squared_z,
                                      z_squared,
                                  })
    );
    EXPECT_EQ(
        graded_lexicographic.GetTerms(), (std::vector<TermType>{
                                             x_squared_z_squared,
                                             x_y_squared_z,
                                             negative_x_cubed,
                                             z_squared,
                                         })
    );
    EXPECT_EQ(
        graded_reverse.GetTerms(), (std::vector<TermType>{
                                       x_y_squared_z,
                                       x_squared_z_squared,
                                       negative_x_cubed,
                                       z_squared,
                                   })
    );

    EXPECT_EQ(lexicographic.GetLeadingCoefficient(), Coefficient{-5});
    EXPECT_EQ(lexicographic.GetLeadingMonomial(), negative_x_cubed.GetMonomial());
    EXPECT_EQ(lexicographic.GetLeadingTerm(), negative_x_cubed);

    EXPECT_EQ(graded_lexicographic.GetLeadingCoefficient(), Coefficient{7});
    EXPECT_EQ(graded_lexicographic.GetLeadingMonomial(), x_squared_z_squared.GetMonomial());
    EXPECT_EQ(graded_lexicographic.GetLeadingTerm(), x_squared_z_squared);

    EXPECT_EQ(graded_reverse.GetLeadingCoefficient(), Coefficient{4});
    EXPECT_EQ(graded_reverse.GetLeadingMonomial(), x_y_squared_z.GetMonomial());
    EXPECT_EQ(graded_reverse.GetLeadingTerm(), x_y_squared_z);

    std::ostringstream lexicographic_output;
    std::ostringstream graded_lexicographic_output;
    std::ostringstream graded_reverse_output;
    lexicographic_output << lexicographic;
    graded_lexicographic_output << graded_lexicographic;
    graded_reverse_output << graded_reverse;

    EXPECT_EQ(lexicographic_output.str(), "-5*x_0^3 + 7*x_0^2*x_2^2 + 4*x_0*x_1^2*x_2 + 4*x_2^2");
    EXPECT_EQ(
        graded_lexicographic_output.str(), "7*x_0^2*x_2^2 + 4*x_0*x_1^2*x_2 - 5*x_0^3 + 4*x_2^2"
    );
    EXPECT_EQ(graded_reverse_output.str(), "4*x_0*x_1^2*x_2 + 7*x_0^2*x_2^2 - 5*x_0^3 + 4*x_2^2");
}

TEST(PolynomialTest, TextbookVectorIsASyzygyOfTheLeadingTerms) {
    const PolynomialType f1{TermType{Coefficient{1}, Monomial{1}}};
    const PolynomialType f2{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
    };
    const PolynomialType f3{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
    };
    const PolynomialType first_component{
        TermType{Coefficient{-1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };
    const PolynomialType second_component{TermType{Coefficient{1}}};
    const PolynomialType third_component{TermType{Coefficient{-1}, Monomial{1}}};

    EXPECT_EQ(f1.GetLeadingTerm(), (TermType{Coefficient{1}, Monomial{1}}));
    EXPECT_EQ(f2.GetLeadingTerm(), (TermType{Coefficient{1}, Monomial{2}}));
    EXPECT_EQ(f3.GetLeadingTerm(), (TermType{Coefficient{1}, Monomial{0, 1}}));

    const PolynomialType leading_term_combination = first_component * f1.GetLeadingTerm() +
                                                    second_component * f2.GetLeadingTerm() +
                                                    third_component * f3.GetLeadingTerm();
    EXPECT_TRUE(leading_term_combination.IsZero());

    const PolynomialType full_polynomial_combination =
        first_component * f1 + second_component * f2 + third_component * f3;
    const PolynomialType expected_full_combination{
        TermType{Coefficient{-1}, Monomial{1, 0, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
    };
    EXPECT_EQ(full_polynomial_combination, expected_full_combination);
}

TEST(PolynomialTest, OneMonicAndZeroPredicatesDistinguishTheirCases) {
    const PolynomialType zero;
    const PolynomialType one{TermType{Coefficient{1}}};
    const PolynomialType other_constant{TermType{Coefficient{2}}};
    const PolynomialType monic_nonconstant{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}},
    };

    EXPECT_TRUE(zero.IsZero());
    EXPECT_FALSE(zero.IsOne());
    EXPECT_FALSE(zero.IsMonic());

    EXPECT_FALSE(one.IsZero());
    EXPECT_TRUE(one.IsOne());
    EXPECT_TRUE(one.IsMonic());

    EXPECT_FALSE(other_constant.IsOne());
    EXPECT_FALSE(other_constant.IsMonic());
    EXPECT_FALSE(monic_nonconstant.IsOne());
    EXPECT_TRUE(monic_nonconstant.IsMonic());
}

TEST(PolynomialTest, MakeMonicDividesEveryCoefficientByTheLeadingCoefficient) {
    PolynomialType polynomial{
        TermType{Coefficient{-2}, Monomial{2}},
        TermType{Coefficient{3}, Monomial{1, 1}},
        TermType{Coefficient{-1}},
    };
    const PolynomialType expected{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-3, 2}, Monomial{1, 1}},
        TermType{Coefficient{1, 2}},
    };

    EXPECT_EQ(&(polynomial.MakeMonic()), &polynomial);
    EXPECT_EQ(polynomial, expected);
    EXPECT_TRUE(polynomial.IsMonic());
}

TEST(PolynomialTest, MakeMonicIsANoOpForZeroAndMonicPolynomials) {
    PolynomialType zero;
    PolynomialType monic{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{-2}},
    };
    const auto expected_monic = monic;

    EXPECT_EQ(&(zero.MakeMonic()), &zero);
    EXPECT_TRUE(zero.IsZero());
    EXPECT_EQ(&(monic.MakeMonic()), &monic);
    EXPECT_EQ(monic, expected_monic);
}

TEST(PolynomialTest, UnaryMinusNegatesAllCoefficientsWithoutChangingTheOperand) {
    const PolynomialType polynomial{
        TermType{Coefficient{2, 3}, Monomial{2}},
        TermType{Coefficient{-4}, Monomial{0, 1}},
    };
    const PolynomialType expected{
        TermType{Coefficient{-2, 3}, Monomial{2}},
        TermType{Coefficient{4}, Monomial{0, 1}},
    };

    EXPECT_EQ(-polynomial, expected);
    EXPECT_EQ(-(-polynomial), polynomial);
    EXPECT_EQ(-PolynomialType{}, PolynomialType{});
}

TEST(PolynomialTest, AdditionSupportsCoefficientsTermsAndPolynomialsOnEitherSide) {
    const PolynomialType polynomial{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{-3}},
    };
    const TermType term{Coefficient{-2}, Monomial{1}};
    const PolynomialType other{
        TermType{Coefficient{5}, Monomial{0, 1}},
        TermType{Coefficient{3}},
    };
    const PolynomialType expected_sum{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{5}, Monomial{0, 1}},
    };

    EXPECT_EQ(polynomial + Coefficient{3}, (PolynomialType{TermType{Coefficient{2}, Monomial{1}}}));
    EXPECT_EQ(Coefficient{3} + polynomial, polynomial + Coefficient{3});
    EXPECT_EQ(polynomial + term, (PolynomialType{TermType{Coefficient{-3}}}));
    EXPECT_EQ(term + polynomial, polynomial + term);
    EXPECT_EQ(polynomial + other, expected_sum);
    EXPECT_EQ(
        polynomial, (PolynomialType{
                        TermType{Coefficient{2}, Monomial{1}},
                        TermType{Coefficient{-3}},
                    })
    );
}

TEST(PolynomialTest, AdditionAssignmentHandlesSelfAdditionAndCancellation) {
    PolynomialType polynomial{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{-3}},
    };
    const PolynomialType doubled{
        TermType{Coefficient{4}, Monomial{1}},
        TermType{Coefficient{-6}},
    };

    EXPECT_EQ(&(polynomial += polynomial), &polynomial);
    EXPECT_EQ(polynomial, doubled);
    EXPECT_EQ(&(polynomial += Coefficient{6}), &polynomial);
    EXPECT_EQ(&(polynomial += TermType{Coefficient{-4}, Monomial{1}}), &polynomial);
    EXPECT_TRUE(polynomial.IsZero());
}

TEST(PolynomialTest, SubtractionSupportsCoefficientsTermsAndPolynomialsOnEitherSide) {
    const PolynomialType polynomial{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{3}},
    };
    const TermType term{Coefficient{1}, Monomial{0, 1}};
    const PolynomialType other{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{-1}},
    };

    EXPECT_EQ(polynomial - Coefficient{3}, (PolynomialType{TermType{Coefficient{2}, Monomial{1}}}));
    EXPECT_EQ(
        Coefficient{3} - polynomial, (PolynomialType{TermType{Coefficient{-2}, Monomial{1}}})
    );
    EXPECT_EQ(
        polynomial - term, (PolynomialType{
                               TermType{Coefficient{2}, Monomial{1}},
                               TermType{Coefficient{-1}, Monomial{0, 1}},
                               TermType{Coefficient{3}},
                           })
    );
    EXPECT_EQ(term - polynomial, -(polynomial - term));
    EXPECT_EQ(polynomial - other, (PolynomialType{TermType{Coefficient{4}}}));
}

TEST(PolynomialTest, SubtractionAssignmentHandlesSelfSubtraction) {
    PolynomialType polynomial{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{3}},
    };

    EXPECT_EQ(&(polynomial -= Coefficient{1}), &polynomial);
    EXPECT_EQ(&(polynomial -= TermType{Coefficient{1}, Monomial{1}}), &polynomial);
    EXPECT_EQ(&(polynomial -= polynomial), &polynomial);
    EXPECT_TRUE(polynomial.IsZero());
}

TEST(PolynomialTest, ScalarMultiplicationHandlesAliasesIdentityAndZero) {
    PolynomialType polynomial{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{-3}},
    };
    const PolynomialType original = polynomial;
    const Coefficient& aliased_factor = polynomial.GetLeadingCoefficient();
    const PolynomialType expected{
        TermType{Coefficient{4}, Monomial{1}},
        TermType{Coefficient{-6}},
    };

    EXPECT_EQ(polynomial * Coefficient{1}, original);
    EXPECT_EQ(Coefficient{2} * polynomial, expected);
    EXPECT_EQ(&(polynomial *= aliased_factor), &polynomial);
    EXPECT_EQ(polynomial, expected);
    EXPECT_TRUE((polynomial * Coefficient{0}).IsZero());
}

TEST(PolynomialTest, TermMultiplicationScalesAndShiftsEveryTerm) {
    const PolynomialType polynomial{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{1}},
    };
    const TermType factor{Coefficient{-3, 2}, Monomial{0, 2}};
    const PolynomialType expected{
        TermType{Coefficient{-3}, Monomial{1, 2}},
        TermType{Coefficient{-3, 2}, Monomial{0, 2}},
    };

    EXPECT_EQ(polynomial * factor, expected);
    EXPECT_EQ(factor * polynomial, expected);
    EXPECT_TRUE((polynomial * TermType{}).IsZero());

    PolynomialType assigned = polynomial;
    EXPECT_EQ(&(assigned *= factor), &assigned);
    EXPECT_EQ(assigned, expected);
}

TEST(PolynomialTest, PolynomialMultiplicationDistributesAndCombinesLikeTerms) {
    const PolynomialType sum{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };
    const PolynomialType difference{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
    };
    const PolynomialType expected{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-1}, Monomial{0, 2}},
    };

    EXPECT_EQ(sum * difference, expected);
    EXPECT_EQ(difference * sum, expected);
    EXPECT_EQ(
        sum, (PolynomialType{
                 TermType{Coefficient{1}, Monomial{1}},
                 TermType{Coefficient{1}, Monomial{0, 1}},
             })
    );
}

TEST(PolynomialTest, PolynomialMultiplicationHandlesZeroOneAndSelfAliasing) {
    PolynomialType polynomial{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };
    const PolynomialType one{TermType{Coefficient{1}}};
    const PolynomialType zero;
    const PolynomialType square{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{2}, Monomial{1, 1}},
        TermType{Coefficient{1}, Monomial{0, 2}},
    };

    EXPECT_EQ(polynomial * one, polynomial);
    EXPECT_EQ(one * polynomial, polynomial);
    EXPECT_TRUE((polynomial * zero).IsZero());
    EXPECT_TRUE((zero * polynomial).IsZero());
    EXPECT_EQ(&(polynomial *= polynomial), &polynomial);
    EXPECT_EQ(polynomial, square);
}

TEST(PolynomialTest, EqualityUsesTheCanonicalTermCollection) {
    const PolynomialType lhs{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{2}, Monomial{0, 1}},
    };
    const PolynomialType reordered_and_split{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{0}, Monomial{8}},
    };

    EXPECT_EQ(lhs, reordered_and_split);
    EXPECT_NE(lhs, (PolynomialType{TermType{Coefficient{1}, Monomial{1}}}));
}

TEST(PolynomialTest, StreamOutputFormatsSignsConstantsAndSpecialCoefficients) {
    const PolynomialType polynomial{
        TermType{Coefficient{-1}, Monomial{2}},
        TermType{Coefficient{3, 2}, Monomial{1, 1}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
        TermType{Coefficient{2}},
    };
    std::ostringstream output;

    output << polynomial;

    EXPECT_EQ(output.str(), "-x_0^2 + 3/2*x_0*x_1 - x_1 + 2");
}

}  // namespace
}  // namespace groebner_basis
