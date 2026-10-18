#include <gtest/gtest.h>

#include <sstream>

#include "groebner_basis/fields/rational.h"

namespace groebner_basis::fields {
namespace {

TEST(RationalTest, DefaultAndIntegerConstructionUseCanonicalRepresentation) {
    const Rational zero;
    const Rational integer{-7};

    EXPECT_TRUE(zero.IsZero());
    EXPECT_EQ(zero.GetNumerator(), Rational::Integer{0});
    EXPECT_EQ(zero.GetDenominator(), Rational::Integer{1});
    EXPECT_FALSE(integer.IsZero());
    EXPECT_EQ(integer.GetNumerator(), Rational::Integer{-7});
    EXPECT_EQ(integer.GetDenominator(), Rational::Integer{1});
}

TEST(RationalTest, FractionConstructionReducesAndNormalizesTheSign) {
    EXPECT_EQ(Rational(12, 18), Rational(2, 3));
    EXPECT_EQ(Rational(-12, 18), Rational(-2, 3));
    EXPECT_EQ(Rational(12, -18), Rational(-2, 3));
    EXPECT_EQ(Rational(-12, -18), Rational(2, 3));
    EXPECT_EQ(Rational(0, -25), Rational(0));
}

#ifndef NDEBUG
TEST(RationalDeathTest, ConstructionRejectsAZeroDenominator) {
    EXPECT_DEATH(
        {
            const Rational invalid(1, 0);
            (void)invalid;
        },
        "denominator must not be zero"
    );
}
#endif

TEST(RationalTest, UnaryMinusNegatesOnlyTheNumerator) {
    const Rational number{3, 4};

    EXPECT_EQ(-number, Rational(-3, 4));
    EXPECT_EQ(-Rational(-3, 4), number);
    EXPECT_EQ(-Rational(0), Rational(0));
    EXPECT_EQ(number, Rational(3, 4));
}

TEST(RationalTest, AdditionAndSubtractionHandleDifferentDenominatorsAndCancellation) {
    const Rational lhs{5, 6};
    const Rational rhs{-7, 15};

    EXPECT_EQ(lhs + rhs, Rational(11, 30));
    EXPECT_EQ(lhs - rhs, Rational(13, 10));
    EXPECT_EQ(lhs + (-lhs), Rational(0));
    EXPECT_EQ(Rational(1, 6) + Rational(5, 6), Rational(1));
    EXPECT_EQ(lhs, Rational(5, 6));
    EXPECT_EQ(rhs, Rational(-7, 15));
}

TEST(RationalTest, AdditionAndSubtractionAssignmentReturnAndModifyTheLeftOperand) {
    Rational value{1, 2};

    Rational* result = &(value += Rational{1, 3});
    EXPECT_EQ(result, &value);
    EXPECT_EQ(value, Rational(5, 6));

    result = &(value -= Rational{7, 6});
    EXPECT_EQ(result, &value);
    EXPECT_EQ(value, Rational(-1, 3));
}

TEST(RationalTest, MultiplicationReducesProductsAndHandlesSignsAndZero) {
    const Rational lhs{-14, 15};
    const Rational rhs{25, 21};

    EXPECT_EQ(lhs * rhs, Rational(-10, 9));
    EXPECT_EQ(lhs * Rational(0), Rational(0));
    EXPECT_EQ(lhs * Rational(1), lhs);
    EXPECT_EQ(lhs, Rational(-14, 15));
    EXPECT_EQ(rhs, Rational(25, 21));
}

TEST(RationalTest, MultiplicationAssignmentReturnsAndModifiesTheLeftOperand) {
    Rational value{4, 9};

    Rational* result = &(value *= Rational{-3, 2});

    EXPECT_EQ(result, &value);
    EXPECT_EQ(value, Rational(-2, 3));
}

TEST(RationalTest, DivisionMultipliesByTheReciprocalAndNormalizesSigns) {
    const Rational lhs{-5, 6};
    const Rational rhs{10, 9};

    EXPECT_EQ(lhs / rhs, Rational(-3, 4));
    EXPECT_EQ(lhs / Rational(-2), Rational(5, 12));
    EXPECT_EQ(Rational(0) / rhs, Rational(0));
    EXPECT_EQ(lhs / Rational(1), lhs);
}

TEST(RationalTest, DivisionAssignmentReturnsAndModifiesTheLeftOperand) {
    Rational value{7, 12};

    Rational* result = &(value /= Rational{-14, 9});

    EXPECT_EQ(result, &value);
    EXPECT_EQ(value, Rational(-3, 8));
}

#ifndef NDEBUG
TEST(RationalDeathTest, DivisionRejectsZero) {
    EXPECT_DEATH(
        {
            Rational value{1};
            value /= Rational{0};
        },
        "divisor must not be zero"
    );
}
#endif

TEST(RationalTest, EqualityAndOrderingUseTheNumericValue) {
    EXPECT_EQ(Rational(2, 4), Rational(1, 2));
    EXPECT_NE(Rational(1, 2), Rational(2, 3));

    EXPECT_LT(Rational(-3, 4), Rational(-2, 3));
    EXPECT_LT(Rational(-1, 2), Rational(0));
    EXPECT_LT(Rational(2, 5), Rational(1, 2));
    EXPECT_LE(Rational(1, 2), Rational(2, 4));
    EXPECT_GT(Rational(7, 6), Rational(1));
    EXPECT_GE(Rational(2, 4), Rational(1, 2));
}

TEST(RationalTest, ArithmeticSatisfiesBasicFieldIdentities) {
    const Rational a{2, 3};
    const Rational b{-5, 7};
    const Rational c{11, 13};

    EXPECT_EQ((a + b) + c, a + (b + c));
    EXPECT_EQ(a * (b + c), a * b + a * c);
    EXPECT_EQ(a - a, Rational(0));
    EXPECT_EQ(a / a, Rational(1));
}

TEST(RationalTest, StreamOutputUsesIntegerOrReducedFractionNotation) {
    std::ostringstream output;

    output << Rational{6, 3} << ' ' << Rational{-6, 8} << ' ' << Rational{0};

    EXPECT_EQ(output.str(), "2 -3/4 0");
}

}  // namespace
}  // namespace groebner_basis::fields
