#include <gtest/gtest.h>

#include <sstream>
#include <utility>
#include <vector>

#include "groebner_basis/monomial.h"

namespace groebner_basis {
namespace {

TEST(MonomialTest, DefaultConstructorCreatesMultiplicativeIdentity) {
    const Monomial monomial;

    EXPECT_TRUE(monomial.IsOne());
    EXPECT_EQ(monomial.GetTotalDegree(), Monomial::Degree{0});
    EXPECT_EQ(monomial.GetVariableIndexUpperBound(), Monomial::VariableIndex{0});
    EXPECT_EQ(monomial.GetDegreeOf(0), Monomial::Degree{0});
    EXPECT_EQ(monomial.GetDegreeOf(100), Monomial::Degree{0});
}

TEST(MonomialTest, InitializerListConstructorStoresDegreesAndTotalDegree) {
    const Monomial monomial{1, 0, 2, 4};

    EXPECT_FALSE(monomial.IsOne());
    EXPECT_EQ(monomial.GetVariableIndexUpperBound(), Monomial::VariableIndex{4});
    EXPECT_EQ(monomial.GetDegreeOf(0), Monomial::Degree{1});
    EXPECT_EQ(monomial.GetDegreeOf(1), Monomial::Degree{0});
    EXPECT_EQ(monomial.GetDegreeOf(2), Monomial::Degree{2});
    EXPECT_EQ(monomial.GetDegreeOf(3), Monomial::Degree{4});
    EXPECT_EQ(monomial.GetDegreeOf(4), Monomial::Degree{0});
    EXPECT_EQ(monomial.GetTotalDegree(), Monomial::Degree{7});
}

TEST(MonomialTest, VectorConstructorCopiesAnLvalueVector) {
    const std::vector<Monomial::Degree> degrees{2, 0, 3};
    const Monomial monomial{degrees};

    EXPECT_EQ(degrees, (std::vector<Monomial::Degree>{2, 0, 3}));
    EXPECT_EQ(monomial.GetDegreeOf(0), Monomial::Degree{2});
    EXPECT_EQ(monomial.GetDegreeOf(1), Monomial::Degree{0});
    EXPECT_EQ(monomial.GetDegreeOf(2), Monomial::Degree{3});
    EXPECT_EQ(monomial.GetTotalDegree(), Monomial::Degree{5});
}

TEST(MonomialTest, ConstructorRemovesOnlyTrailingZeroDegrees) {
    const Monomial with_trailing_zeros{1, 0, 2, 0, 0};
    const Monomial expected{1, 0, 2};

    EXPECT_EQ(with_trailing_zeros, expected);
    EXPECT_EQ(with_trailing_zeros.GetVariableIndexUpperBound(), Monomial::VariableIndex{3});
    EXPECT_EQ(with_trailing_zeros.GetTotalDegree(), Monomial::Degree{3});

    const Monomial with_internal_zero{0, 2, 0};
    EXPECT_EQ(with_internal_zero.GetVariableIndexUpperBound(), Monomial::VariableIndex{2});
    EXPECT_EQ(with_internal_zero.GetDegreeOf(0), Monomial::Degree{0});
    EXPECT_EQ(with_internal_zero.GetDegreeOf(1), Monomial::Degree{2});
}

TEST(MonomialTest, AllZeroDegreesCreateMultiplicativeIdentity) {
    const Monomial all_zeros{0, 0, 0};

    EXPECT_TRUE(all_zeros.IsOne());
    EXPECT_EQ(all_zeros, Monomial{});
    EXPECT_EQ(all_zeros.GetVariableIndexUpperBound(), Monomial::VariableIndex{0});
    EXPECT_EQ(all_zeros.GetTotalDegree(), Monomial::Degree{0});
}

TEST(MonomialTest, EqualityUsesCanonicalDegreeRepresentation) {
    const Monomial lhs{1, 0, 2, 0};
    const Monomial equal{1, 0, 2};
    const Monomial different_degree{1, 1, 2};
    const Monomial different_variable{1, 0, 0, 2};

    EXPECT_EQ(lhs, equal);
    EXPECT_NE(lhs, different_degree);
    EXPECT_NE(lhs, different_variable);
}

TEST(MonomialTest, IdentityDividesEveryMonomial) {
    const Monomial identity;
    const Monomial monomial{4, 0, 2};

    EXPECT_TRUE(identity.Divides(identity));
    EXPECT_TRUE(identity.Divides(monomial));
    EXPECT_FALSE(monomial.Divides(identity));
}

TEST(MonomialTest, DividesComparesEveryVariableDegree) {
    const Monomial divisor{1, 0, 2};
    const Monomial dividend{3, 4, 2, 1};
    const Monomial excessive_first_degree{4, 0, 2};
    const Monomial excessive_last_degree{1, 0, 3};
    const Monomial extra_variable{0, 0, 0, 0, 1};

    EXPECT_TRUE(divisor.Divides(divisor));
    EXPECT_TRUE(divisor.Divides(dividend));
    EXPECT_FALSE(excessive_first_degree.Divides(dividend));
    EXPECT_FALSE(excessive_last_degree.Divides(dividend));
    EXPECT_FALSE(extra_variable.Divides(dividend));
}

TEST(MonomialTest, MultiplicationAddsDegreesAndPreservesOperands) {
    const Monomial lhs{1, 0, 2};
    const Monomial rhs{0, 3, 1, 4};
    const Monomial expected{1, 3, 3, 4};

    const Monomial product = lhs * rhs;

    EXPECT_EQ(product, expected);
    EXPECT_EQ(product.GetTotalDegree(), Monomial::Degree{11});
    EXPECT_EQ(lhs, (Monomial{1, 0, 2}));
    EXPECT_EQ(rhs, (Monomial{0, 3, 1, 4}));
}

TEST(MonomialTest, MultiplicationAssignmentSupportsIdentityAndSelfMultiplication) {
    const Monomial identity;
    Monomial monomial{1, 0, 2};

    monomial *= identity;
    EXPECT_EQ(monomial, (Monomial{1, 0, 2}));
    EXPECT_EQ(identity * monomial, monomial);

    monomial *= monomial;
    EXPECT_EQ(monomial, (Monomial{2, 0, 4}));
    EXPECT_EQ(monomial.GetTotalDegree(), Monomial::Degree{6});
}

TEST(MonomialTest, DivisionSubtractsDegreesAndPreservesOperands) {
    const Monomial dividend{3, 4, 2};
    const Monomial divisor{1, 4, 1};
    const Monomial expected{2, 0, 1};

    const Monomial quotient = dividend / divisor;

    EXPECT_EQ(quotient, expected);
    EXPECT_EQ(quotient.GetTotalDegree(), Monomial::Degree{3});
    EXPECT_EQ(dividend, (Monomial{3, 4, 2}));
    EXPECT_EQ(divisor, (Monomial{1, 4, 1}));
}

TEST(MonomialTest, DivisionAssignmentCanonicalizesTheResult) {
    Monomial quotient{1, 3};
    const Monomial divisor{0, 3};

    quotient /= divisor;

    EXPECT_EQ(quotient, (Monomial{1}));
    EXPECT_EQ(quotient.GetVariableIndexUpperBound(), Monomial::VariableIndex{1});
    EXPECT_EQ(quotient.GetTotalDegree(), Monomial::Degree{1});
}

TEST(MonomialTest, DivisionByIdentityAndSelfDivisionProduceExpectedResults) {
    const Monomial identity;
    Monomial monomial{2, 1, 3};

    EXPECT_EQ(monomial / identity, monomial);

    monomial /= monomial;
    EXPECT_TRUE(monomial.IsOne());
    EXPECT_EQ(monomial.GetTotalDegree(), Monomial::Degree{0});
}

#ifndef NDEBUG
TEST(MonomialDeathTest, DivisionRejectsANondivisor) {
    EXPECT_DEATH(
        {
            Monomial dividend(std::vector<Monomial::Degree>{1, 1});
            dividend /= Monomial{2};
        },
        "divisor divides dividend");
}
#endif

TEST(MonomialTest, CopyOperationsCreateIndependentValues) {
    const Monomial original{1, 2, 3};
    Monomial copy{original};

    copy *= Monomial{1};

    EXPECT_EQ(original, (Monomial{1, 2, 3}));
    EXPECT_EQ(copy, (Monomial{2, 2, 3}));

    Monomial assigned;
    assigned = original;
    assigned /= Monomial{0, 1, 1};

    EXPECT_EQ(original, (Monomial{1, 2, 3}));
    EXPECT_EQ(assigned, (Monomial{1, 1, 2}));
}

TEST(MonomialTest, MoveOperationsPreserveTheTransferredValue) {
    Monomial move_source{1, 2, 3};
    Monomial move_constructed{std::move(move_source)};

    EXPECT_EQ(move_constructed, (Monomial{1, 2, 3}));
    EXPECT_EQ(move_constructed.GetTotalDegree(), Monomial::Degree{6});

    Monomial assignment_source{0, 4, 1};
    Monomial move_assigned{7};
    move_assigned = std::move(assignment_source);

    EXPECT_EQ(move_assigned, (Monomial{0, 4, 1}));
    EXPECT_EQ(move_assigned.GetTotalDegree(), Monomial::Degree{5});
}

TEST(MonomialTest, SwapExchangesDegreesAndCachedTotalDegree) {
    Monomial lhs{1, 2};
    Monomial rhs{0, 0, 4};

    swap(lhs, rhs);

    EXPECT_EQ(lhs, (Monomial{0, 0, 4}));
    EXPECT_EQ(lhs.GetTotalDegree(), Monomial::Degree{4});
    EXPECT_EQ(rhs, (Monomial{1, 2}));
    EXPECT_EQ(rhs.GetTotalDegree(), Monomial::Degree{3});
}

TEST(MonomialTest, StreamOutputFormatsIdentityVariablesAndPowers) {
    std::ostringstream output;

    output << Monomial{};
    EXPECT_EQ(output.str(), "1");

    output.str("");
    output.clear();
    output << Monomial{1};
    EXPECT_EQ(output.str(), "x_0");

    output.str("");
    output.clear();
    output << Monomial{0, 3};
    EXPECT_EQ(output.str(), "x_1^3");

    output.str("");
    output.clear();
    output << Monomial{1, 0, 2};
    EXPECT_EQ(output.str(), "x_0*x_2^2");
}

TEST(MonomialTest, LeastCommonMultipleUsesMaximumDegreeForEachVariable) {
    const Monomial lhs{1, 0, 2};
    const Monomial rhs{0, 3, 1, 4};
    const Monomial expected{1, 3, 2, 4};

    const Monomial least_common_multiple = ComputeLeastCommonMultiple(lhs, rhs);

    EXPECT_EQ(least_common_multiple, expected);
    EXPECT_EQ(least_common_multiple.GetTotalDegree(), Monomial::Degree{10});
    EXPECT_TRUE(lhs.Divides(least_common_multiple));
    EXPECT_TRUE(rhs.Divides(least_common_multiple));
    EXPECT_EQ(ComputeLeastCommonMultiple(rhs, lhs), least_common_multiple);
}

TEST(MonomialTest, LeastCommonMultipleHandlesIdentityAndEqualOperands) {
    const Monomial identity;
    const Monomial monomial{2, 0, 5};

    EXPECT_EQ(ComputeLeastCommonMultiple(identity, monomial), monomial);
    EXPECT_EQ(ComputeLeastCommonMultiple(monomial, identity), monomial);
    EXPECT_EQ(ComputeLeastCommonMultiple(monomial, monomial), monomial);
    EXPECT_EQ(ComputeLeastCommonMultiple(identity, identity), identity);
}

}  // namespace
}  // namespace groebner_basis
