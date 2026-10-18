#include <gtest/gtest.h>

#include <vector>

#include "groebner_basis/fields/rational.h"
#include "groebner_basis/polynomial_ideal.h"

namespace groebner_basis {
namespace {

using Coefficient = fields::Rational;
using TermType = Term<Coefficient>;
using PolynomialType = Polynomial<Coefficient, LexicographicOrder>;
using Ideal = PolynomialIdeal<Coefficient, LexicographicOrder>;

const PolynomialType kX{TermType{Coefficient{1}, Monomial{1}}};
const PolynomialType kY{TermType{Coefficient{1}, Monomial{0, 1}}};
const PolynomialType kOne{TermType{Coefficient{1}}};

TEST(PolynomialIdealTest, DefaultConstructionCreatesTheZeroIdeal) {
    const Ideal ideal;

    EXPECT_TRUE(ideal.IsZero());
    EXPECT_FALSE(ideal.IsUnit());
    EXPECT_TRUE(ideal.GetReducedGroebnerBasis().empty());
    EXPECT_TRUE(ideal.ContainsPolynomial(PolynomialType{}));
    EXPECT_FALSE(ideal.ContainsPolynomial(kX));
    EXPECT_EQ(ideal.ComputeNormalForm(kX), kX);
}

TEST(PolynomialIdealTest, ConstructionComputesACanonicalReducedGroebnerBasis) {
    const PolynomialType first{
        TermType{Coefficient{2}, Monomial{1}},
        TermType{Coefficient{2}, Monomial{0, 1}},
    };
    const PolynomialType second{
        TermType{Coefficient{-3}, Monomial{1}},
        TermType{Coefficient{3}, Monomial{0, 1}},
    };
    const Ideal ideal{PolynomialType{}, first, second};
    const std::vector<PolynomialType> expected{kX, kY};

    EXPECT_FALSE(ideal.IsZero());
    EXPECT_FALSE(ideal.IsUnit());
    EXPECT_EQ(ideal.GetReducedGroebnerBasis(), expected);
}

TEST(PolynomialIdealTest, VectorConstructionLeavesAnLvalueGeneratorListUnchanged) {
    const std::vector<PolynomialType> generators{kX, kY};
    const auto original = generators;

    const Ideal ideal{generators};

    EXPECT_EQ(generators, original);
    EXPECT_EQ(ideal.GetReducedGroebnerBasis(), original);
}

TEST(PolynomialIdealTest, UnitIdealIsDetectedAndContainsEveryPolynomial) {
    const Ideal unit{
        PolynomialType{TermType{Coefficient{-7}}},
        kX + kY,
    };
    const PolynomialType arbitrary{
        TermType{Coefficient{3, 5}, Monomial{4, 2}},
        TermType{Coefficient{-9}, Monomial{0, 3}},
        TermType{Coefficient{11, 7}},
    };

    ASSERT_EQ(unit.GetReducedGroebnerBasis().size(), 1U);
    EXPECT_TRUE(unit.IsUnit());
    EXPECT_FALSE(unit.IsZero());
    EXPECT_TRUE(unit.GetReducedGroebnerBasis().front().IsOne());
    EXPECT_TRUE(unit.ContainsPolynomial(arbitrary));
    EXPECT_TRUE(unit.ComputeNormalForm(arbitrary).IsZero());
}

TEST(PolynomialIdealTest, MembershipAndNormalFormUseTheReducedBasis) {
    const PolynomialType y_squared{TermType{Coefficient{1}, Monomial{0, 2}}};
    const Ideal ideal{kX, y_squared};
    const PolynomialType member{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{3}, Monomial{0, 2}},
    };
    const PolynomialType nonmember{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{1}, Monomial{0, 3}},
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{2}},
    };
    const PolynomialType expected_normal_form{
        TermType{Coefficient{1}, Monomial{0, 1}},
        TermType{Coefficient{2}},
    };

    EXPECT_TRUE(ideal.ContainsPolynomial(member));
    EXPECT_FALSE(ideal.ContainsPolynomial(nonmember));
    const auto normal_form = ideal.ComputeNormalForm(nonmember);
    EXPECT_EQ(normal_form, expected_normal_form);
    EXPECT_EQ(ideal.ComputeNormalForm(normal_form), normal_form);
}

TEST(PolynomialIdealTest, ExplicitDecompositionProvesFourVariableMembership) {
    const PolynomialType first_generator{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{0, 1}},
    };
    const PolynomialType second_generator{
        TermType{Coefficient{1}, Monomial{0, 0, 1}},
        TermType{Coefficient{1}, Monomial{0, 0, 0, 1}},
    };
    const PolynomialType first_multiplier{TermType{Coefficient{1}, Monomial{0, 0, 2, 1}}};
    const PolynomialType second_multiplier{TermType{Coefficient{1}, Monomial{1, 0, 1, 1}}};
    const PolynomialType polynomial{
        TermType{Coefficient{1}, Monomial{0, 1, 2, 1}},
        TermType{Coefficient{-1}, Monomial{1, 0, 1, 2}},
    };
    const TermType first_reduction_quotient{
        Coefficient{-1},
        Monomial{0, 0, 1, 2},
    };
    const PolynomialType after_first_reduction{
        TermType{Coefficient{1}, Monomial{0, 1, 2, 1}},
        TermType{Coefficient{1}, Monomial{0, 1, 1, 2}},
    };
    const PolynomialType two_reduction_quotient{
        TermType{Coefficient{1}, Monomial{0, 1, 1, 1}},
        TermType{Coefficient{-1}, Monomial{0, 1, 0, 2}},
    };
    const PolynomialType after_two_more_reductions{
        TermType{Coefficient{1}, Monomial{0, 1, 1, 2}},
        TermType{Coefficient{1}, Monomial{0, 1, 0, 3}},
    };
    const TermType final_reduction_quotient{
        Coefficient{1},
        Monomial{0, 1, 0, 2},
    };
    const Ideal ideal{first_generator, second_generator};

    EXPECT_EQ(
        polynomial, first_multiplier * first_generator - second_multiplier * second_generator
    );
    EXPECT_EQ(polynomial - first_reduction_quotient * first_generator, after_first_reduction);
    EXPECT_EQ(
        after_first_reduction - two_reduction_quotient * second_generator, after_two_more_reductions
    );
    EXPECT_TRUE((after_two_more_reductions - final_reduction_quotient * second_generator).IsZero());

    const auto direct_second_factor =
        after_first_reduction.GetLeadingTerm() / second_generator.GetLeadingTerm();
    EXPECT_EQ(direct_second_factor, (TermType{Coefficient{1}, Monomial{0, 1, 1, 1}}));
    EXPECT_TRUE((after_first_reduction - direct_second_factor * second_generator).IsZero());
    EXPECT_TRUE(
        algorithms::ReducePolynomial(
            polynomial, std::vector<PolynomialType>{first_generator, second_generator}
        )
            .IsZero()
    );
    EXPECT_TRUE(ideal.ContainsPolynomial(polynomial));
    EXPECT_TRUE(ideal.ComputeNormalForm(polynomial).IsZero());
}

TEST(PolynomialIdealTest, NonzeroNormalFormProvesTwoVariableNonmembership) {
    const PolynomialType x_squared{TermType{Coefficient{1}, Monomial{2}}};
    const PolynomialType y{TermType{Coefficient{1}, Monomial{0, 1}}};
    const PolynomialType x{TermType{Coefficient{1}, Monomial{1}}};
    const PolynomialType polynomial{
        TermType{Coefficient{1}, Monomial{1}},
        TermType{Coefficient{1}, Monomial{1, 2}},
        TermType{Coefficient{3}, Monomial{1, 3}},
    };
    const PolynomialType y_multiplier{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{3}, Monomial{1, 2}},
    };
    const Ideal ideal{x_squared, y};

    EXPECT_EQ(polynomial, x + y * y_multiplier);
    EXPECT_FALSE(x_squared.GetLeadingMonomial().Divides(x.GetLeadingMonomial()));
    EXPECT_FALSE(y.GetLeadingMonomial().Divides(x.GetLeadingMonomial()));
    EXPECT_FALSE(ideal.ContainsPolynomial(polynomial));
    EXPECT_EQ(ideal.ComputeNormalForm(polynomial), x);
}

TEST(PolynomialIdealTest, GradedLexBasisDecidesTextbookMembershipExamples) {
    using GradedPolynomial = Polynomial<Coefficient, GradedLexOrder>;
    using GradedIdeal = PolynomialIdeal<Coefficient, GradedLexOrder>;

    const GradedPolynomial f1{
        TermType{Coefficient{1}, Monomial{1, 0, 1}},
        TermType{Coefficient{-1}, Monomial{0, 2}},
    };
    const GradedPolynomial f2{
        TermType{Coefficient{1}, Monomial{3}},
        TermType{Coefficient{-1}, Monomial{0, 0, 2}},
    };
    const GradedPolynomial f3{
        TermType{Coefficient{1}, Monomial{2, 2}},
        TermType{Coefficient{-1}, Monomial{0, 0, 3}},
    };
    const GradedPolynomial f4{
        TermType{Coefficient{1}, Monomial{1, 4}},
        TermType{Coefficient{-1}, Monomial{0, 0, 4}},
    };
    const GradedPolynomial f5{
        TermType{Coefficient{1}, Monomial{0, 6}},
        TermType{Coefficient{-1}, Monomial{0, 0, 5}},
    };
    const std::vector<GradedPolynomial> expected_reduced_basis{
        f5, f4, f3, f2, f1,
    };
    const GradedPolynomial first_s_polynomial{
        TermType{Coefficient{-1}, Monomial{2, 2}},
        TermType{Coefficient{1}, Monomial{0, 0, 3}},
    };

    EXPECT_EQ(algorithms::ComputeSPolynomial(f1, f2), first_s_polynomial);
    EXPECT_EQ(
        algorithms::ReducePolynomial(first_s_polynomial, std::vector<GradedPolynomial>{f1, f2}),
        first_s_polynomial
    );
    EXPECT_FALSE(f1.GetLeadingMonomial().Divides(first_s_polynomial.GetLeadingMonomial()));
    EXPECT_FALSE(f2.GetLeadingMonomial().Divides(first_s_polynomial.GetLeadingMonomial()));

    const GradedIdeal ideal{f1, f2};
    EXPECT_EQ(ideal.GetReducedGroebnerBasis(), expected_reduced_basis);

    for (const auto& polynomial : expected_reduced_basis) {
        EXPECT_TRUE(polynomial.IsMonic());
        for (const auto& term : polynomial.GetTerms()) {
            if (term.GetMonomial() == polynomial.GetLeadingMonomial()) {
                continue;
            }
            for (const auto& other : expected_reduced_basis) {
                EXPECT_FALSE(other.GetLeadingMonomial().Divides(term.GetMonomial()));
            }
        }
    }

    const GradedPolynomial member{
        TermType{Coefficient{-4}, Monomial{2, 2, 2}},
        TermType{Coefficient{1}, Monomial{0, 6}},
        TermType{Coefficient{3}, Monomial{0, 0, 5}},
    };
    const GradedPolynomial first_quotient{
        TermType{Coefficient{-4}, Monomial{1, 2, 1}},
        TermType{Coefficient{-4}, Monomial{0, 4}},
    };

    EXPECT_EQ(member, first_quotient * f1 + Coefficient{-3} * f5);
    EXPECT_TRUE(ideal.ContainsPolynomial(member));
    EXPECT_TRUE(ideal.ComputeNormalForm(member).IsZero());

    const GradedPolynomial nonmember{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-5}, Monomial{0, 0, 2}},
        TermType{Coefficient{1}, Monomial{1}},
    };

    EXPECT_EQ(nonmember.GetLeadingMonomial(), (Monomial{1, 1}));
    for (const auto& polynomial : expected_reduced_basis) {
        EXPECT_FALSE(polynomial.GetLeadingMonomial().Divides(nonmember.GetLeadingMonomial()));
    }
    EXPECT_FALSE(ideal.ContainsPolynomial(nonmember));
    EXPECT_EQ(ideal.ComputeNormalForm(nonmember), nonmember);
}

TEST(PolynomialIdealTest, EqualityIsIndependentOfGeneratorChoiceAndScaling) {
    const Ideal coordinate_ideal{kX, kY};
    const Ideal alternate_generators{
        PolynomialType{
            TermType{Coefficient{2}, Monomial{1}},
            TermType{Coefficient{2}, Monomial{0, 1}},
        },
        PolynomialType{
            TermType{Coefficient{3}, Monomial{1}},
            TermType{Coefficient{-3}, Monomial{0, 1}},
        },
    };
    const Ideal x_axis{kX};

    EXPECT_EQ(coordinate_ideal, alternate_generators);
    EXPECT_NE(coordinate_ideal, x_axis);
}

TEST(PolynomialIdealTest, AdditionComputesTheIdealGeneratedByBothOperands) {
    const PolynomialType x_squared{TermType{Coefficient{1}, Monomial{2}}};
    const Ideal smaller{x_squared};
    const Ideal larger{kX};
    const Ideal zero;

    EXPECT_EQ(smaller + larger, larger);
    EXPECT_EQ(larger + smaller, larger);
    EXPECT_EQ(smaller + zero, smaller);
    EXPECT_EQ(zero + smaller, smaller);

    Ideal assigned = smaller;
    EXPECT_EQ(&(assigned += larger), &assigned);
    EXPECT_EQ(assigned, larger);
}

TEST(PolynomialIdealTest, SelfAdditionIsANoOp) {
    Ideal ideal{kX, kY};
    const Ideal original = ideal;

    EXPECT_EQ(&(ideal += ideal), &ideal);
    EXPECT_EQ(ideal, original);
}

TEST(PolynomialIdealTest, MultiplicationUsesAllPairwiseGeneratorProducts) {
    const Ideal x_ideal{kX};
    const Ideal y_ideal{kY};
    const Ideal expected{kX * kY};

    EXPECT_EQ(x_ideal * y_ideal, expected);
    EXPECT_EQ(y_ideal * x_ideal, expected);

    Ideal assigned = x_ideal;
    EXPECT_EQ(&(assigned *= y_ideal), &assigned);
    EXPECT_EQ(assigned, expected);
}

TEST(PolynomialIdealTest, MultiplicationHandlesZeroUnitAndSelfAliasing) {
    const Ideal zero;
    const Ideal unit{kOne};
    const Ideal x_ideal{kX};
    const Ideal x_squared_ideal{kX * kX};

    EXPECT_EQ(x_ideal * unit, x_ideal);
    EXPECT_EQ(unit * x_ideal, x_ideal);
    EXPECT_TRUE((x_ideal * zero).IsZero());
    EXPECT_TRUE((zero * x_ideal).IsZero());

    Ideal assigned = x_ideal;
    EXPECT_EQ(&(assigned *= assigned), &assigned);
    EXPECT_EQ(assigned, x_squared_ideal);
}

}  // namespace
}  // namespace groebner_basis
