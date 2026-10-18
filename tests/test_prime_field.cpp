#include <gtest/gtest.h>

#include <cstdint>
#include <sstream>

#include "groebner_basis/fields/prime_field.h"

namespace groebner_basis::fields {
namespace {

using F7 = PrimeField<7>;

static_assert(detail::IsPrime(2));
static_assert(detail::IsPrime(97));
static_assert(!detail::IsPrime(-3));
static_assert(!detail::IsPrime(0));
static_assert(!detail::IsPrime(1));
static_assert(!detail::IsPrime(49));
static_assert(F7::GetModulus() == 7);

TEST(PrimeFieldDetailTest, ModularExponentiationReducesPositiveAndNegativeBases) {
    EXPECT_EQ(detail::ComputeModularExponentiation(3, 0, 7), 1);
    EXPECT_EQ(detail::ComputeModularExponentiation(3, 4, 7), 4);
    EXPECT_EQ(detail::ComputeModularExponentiation(-2, 5, 7), 3);
}

TEST(PrimeFieldTest, DefaultConstructionCreatesZero) {
    const F7 value;

    EXPECT_EQ(value.GetValue(), F7::Integer{0});
    EXPECT_EQ(value, F7{0});
}

TEST(PrimeFieldTest, ConstructionReducesEveryIntegerToItsCanonicalResidue) {
    EXPECT_EQ(F7{3}.GetValue(), F7::Integer{3});
    EXPECT_EQ(F7{10}.GetValue(), F7::Integer{3});
    EXPECT_EQ(F7{-1}.GetValue(), F7::Integer{6});
    EXPECT_EQ(F7{-15}.GetValue(), F7::Integer{6});
    EXPECT_EQ(F7{21}.GetValue(), F7::Integer{0});
}

TEST(PrimeFieldTest, UnaryMinusComputesTheAdditiveInverse) {
    EXPECT_EQ(-F7{3}, F7{4});
    EXPECT_EQ(-F7{-1}, F7{1});
    EXPECT_EQ(-F7{0}, F7{0});
}

TEST(PrimeFieldTest, AdditionWrapsAtTheModulus) {
    EXPECT_EQ(F7{3} + F7{2}, F7{5});
    EXPECT_EQ(F7{5} + F7{4}, F7{2});
    EXPECT_EQ(F7{6} + F7{1}, F7{0});
}

TEST(PrimeFieldTest, SubtractionWrapsBelowZero) {
    EXPECT_EQ(F7{5} - F7{3}, F7{2});
    EXPECT_EQ(F7{2} - F7{5}, F7{4});
    EXPECT_EQ(F7{3} - F7{3}, F7{0});
}

TEST(PrimeFieldTest, MultiplicationReducesTheProduct) {
    EXPECT_EQ(F7{3} * F7{4}, F7{5});
    EXPECT_EQ(F7{-2} * F7{3}, F7{1});
    EXPECT_EQ(F7{0} * F7{6}, F7{0});
}

TEST(PrimeFieldTest, DivisionUsesTheUniqueMultiplicativeInverse) {
    EXPECT_EQ(F7{3} / F7{2}, F7{5});
    EXPECT_EQ(F7{1} / F7{3}, F7{5});
    EXPECT_EQ(F7{0} / F7{6}, F7{0});

    for (std::int64_t value = 1; value < F7::GetModulus(); ++value) {
        const F7 element{value};
        EXPECT_EQ(element / element, F7{1}) << "value = " << value;
        EXPECT_EQ(F7{1} / element * element, F7{1}) << "value = " << value;
    }
}

TEST(PrimeFieldTest, CompoundAssignmentsReturnAndModifyTheLeftOperand) {
    F7 value{5};

    EXPECT_EQ(&(value += F7{4}), &value);
    EXPECT_EQ(value, F7{2});
    EXPECT_EQ(&(value -= F7{5}), &value);
    EXPECT_EQ(value, F7{4});
    EXPECT_EQ(&(value *= F7{3}), &value);
    EXPECT_EQ(value, F7{5});
    EXPECT_EQ(&(value /= F7{2}), &value);
    EXPECT_EQ(value, F7{6});
}

#ifndef NDEBUG
TEST(PrimeFieldDeathTest, DivisionRejectsZero) {
    EXPECT_DEATH(
        {
            F7 value{1};
            value /= F7{0};
        },
        "divisor must not be zero"
    );
}
#endif

TEST(PrimeFieldTest, EqualityAndOrderingCompareCanonicalResidues) {
    EXPECT_EQ(F7{10}, F7{3});
    EXPECT_NE(F7{2}, F7{3});
    EXPECT_LT(F7{2}, F7{3});
    EXPECT_LE(F7{10}, F7{3});
    EXPECT_GT(F7{-1}, F7{3});
    EXPECT_GE(F7{3}, F7{10});
}

TEST(PrimeFieldTest, ArithmeticSatisfiesFieldIdentitiesExhaustively) {
    for (std::int64_t lhs_value = 0; lhs_value < F7::GetModulus(); ++lhs_value) {
        for (std::int64_t rhs_value = 0; rhs_value < F7::GetModulus(); ++rhs_value) {
            const F7 lhs{lhs_value};
            const F7 rhs{rhs_value};

            EXPECT_EQ(lhs + rhs, rhs + lhs);
            EXPECT_EQ(lhs * rhs, rhs * lhs);
            EXPECT_EQ(lhs + F7{0}, lhs);
            EXPECT_EQ(lhs * F7{1}, lhs);
            EXPECT_EQ(lhs + (-lhs), F7{0});
            EXPECT_EQ((lhs + rhs) - rhs, lhs);

            for (std::int64_t third_value = 0; third_value < F7::GetModulus(); ++third_value) {
                const F7 third{third_value};
                EXPECT_EQ(lhs * (rhs + third), lhs * rhs + lhs * third);
            }
        }
    }
}

TEST(PrimeFieldTest, StreamOutputWritesTheCanonicalResidue) {
    std::ostringstream output;

    output << F7{-1} << ' ' << F7{14} << ' ' << F7{3};

    EXPECT_EQ(output.str(), "6 0 3");
}

}  // namespace
}  // namespace groebner_basis::fields
