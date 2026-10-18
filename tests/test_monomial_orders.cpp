#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "groebner_basis/monomial.h"
#include "groebner_basis/monomial_orders.h"

using namespace groebner_basis;

TEST(MonomialOrderTest, LexicographicOrder1) {
    // 581
    Monomial monomial_a = Monomial({3, 2, 1});
    Monomial monomial_b = Monomial({2, 6, 12});

    EXPECT_TRUE(LexicographicOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(LexicographicOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, LexicographicOrder2) {
    // 581
    Monomial monomial_a = Monomial({3, 6});
    Monomial monomial_b = Monomial({3, 4, 1});

    EXPECT_TRUE(LexicographicOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(LexicographicOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, LexicographicOrder3) {
    // 581
    Monomial monomial_a = Monomial({2, 2, 2});
    Monomial monomial_b = Monomial({1, 4, 1});

    EXPECT_TRUE(LexicographicOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(LexicographicOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, LexicographicOrder4) {
    // 653
    Monomial monomial_a = Monomial({1, 2});
    Monomial monomial_b = Monomial({0, 3, 4});

    EXPECT_TRUE(LexicographicOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(LexicographicOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, LexicographicOrder5) {
    // 653
    Monomial monomial_a = Monomial({3, 2, 4});
    Monomial monomial_b = Monomial({3, 2, 1});

    EXPECT_TRUE(LexicographicOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(LexicographicOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, LexicographicOrderUsesTheFirstDifferingVariable) {
    const Monomial monomial_a{1, 2, 4};
    const Monomial monomial_b{1, 1, 5};

    EXPECT_TRUE(LexicographicOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(LexicographicOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, LexicographicTupleExamplesUseTheFirstDifferingCoordinate) {
    const Monomial first{2, 3, 0, 7};
    const Monomial second{4, 0, 0, 0};

    EXPECT_FALSE(LexicographicOrder{}(first, second));
    EXPECT_TRUE(LexicographicOrder{}(second, first));

    const Monomial third{3, 1, 5, 2};
    const Monomial fourth{3, 1, 5, 1};

    EXPECT_TRUE(LexicographicOrder{}(third, fourth));
    EXPECT_FALSE(LexicographicOrder{}(fourth, third));
}

TEST(MonomialOrderTest, LexicographicOrder6) {
    // 581
    Monomial monomial_a = Monomial({3, 2});
    Monomial monomial_b = Monomial({2, 1, 3});

    EXPECT_TRUE(LexicographicOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(LexicographicOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, GradedReverseOrder1) {
    // 581
    Monomial monomial_a = Monomial({2, 6, 12});
    Monomial monomial_b = Monomial({3, 2, 1});

    EXPECT_TRUE(GradedReverseLexOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(GradedReverseLexOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, GradedReverseOrder2) {
    // 581
    Monomial monomial_a = Monomial({1, 4, 1});
    Monomial monomial_b = Monomial({2, 2, 2});

    EXPECT_TRUE(GradedReverseLexOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(GradedReverseLexOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, GradedReverseOrder3) {
    // 653
    Monomial monomial_a = Monomial({4, 7, 1});
    Monomial monomial_b = Monomial({4, 2, 3});

    EXPECT_TRUE(GradedReverseLexOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(GradedReverseLexOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, GradedReverseOrder4) {
    // 653
    Monomial monomial_a = Monomial({1, 5, 2});
    Monomial monomial_b = Monomial({4, 1, 3});

    EXPECT_TRUE(GradedReverseLexOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(GradedReverseLexOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, GradedReverseOrder5) {
    // 581
    Monomial monomial_a = Monomial({3, 2});
    Monomial monomial_b = Monomial({2, 1, 3});

    EXPECT_FALSE(GradedReverseLexOrder{}(monomial_a, monomial_b));
    EXPECT_TRUE(GradedReverseLexOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, GradedLexOrder1) {
    // 653
    Monomial monomial_a = Monomial({1, 2, 3});
    Monomial monomial_b = Monomial({3, 2});

    EXPECT_TRUE(GradedLexOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(GradedLexOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, GradedLexOrder2) {
    // 653
    Monomial monomial_a = Monomial({1, 2, 4});
    Monomial monomial_b = Monomial({1, 1, 5});

    EXPECT_TRUE(GradedLexOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(GradedLexOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, AllOrders1) {
    // 653
    Monomial monomial_a = Monomial({5, 1, 1});
    Monomial monomial_b = Monomial({4, 1, 2});

    EXPECT_TRUE(LexicographicOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(LexicographicOrder{}(monomial_b, monomial_a));

    EXPECT_TRUE(GradedLexOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(GradedLexOrder{}(monomial_b, monomial_a));

    EXPECT_TRUE(GradedReverseLexOrder{}(monomial_a, monomial_b));
    EXPECT_FALSE(GradedReverseLexOrder{}(monomial_b, monomial_a));
}

TEST(MonomialOrderTest, LexicographicOrderFull) {
    std::vector<Monomial> monomials = {
        Monomial({1, 2, 1}), Monomial({0, 0, 2}), Monomial({3}), Monomial({2, 0, 2})
    };

    EXPECT_TRUE(LexicographicOrder{}(monomials[2], monomials[3]));
    EXPECT_TRUE(LexicographicOrder{}(monomials[3], monomials[0]));
    EXPECT_TRUE(LexicographicOrder{}(monomials[0], monomials[1]));

    EXPECT_TRUE(LexicographicOrder{}(monomials[2], monomials[0]));
    EXPECT_TRUE(LexicographicOrder{}(monomials[2], monomials[1]));
    EXPECT_TRUE(LexicographicOrder{}(monomials[3], monomials[1]));

    EXPECT_FALSE(LexicographicOrder{}(monomials[3], monomials[2]));
    EXPECT_FALSE(LexicographicOrder{}(monomials[0], monomials[3]));
    EXPECT_FALSE(LexicographicOrder{}(monomials[1], monomials[0]));

    EXPECT_FALSE(LexicographicOrder{}(monomials[0], monomials[2]));
    EXPECT_FALSE(LexicographicOrder{}(monomials[1], monomials[2]));
    EXPECT_FALSE(LexicographicOrder{}(monomials[1], monomials[3]));
}

TEST(MonomialOrderTest, GradedLexOrderFull) {
    std::vector<Monomial> monomials = {
        Monomial({1, 2, 1}), Monomial({0, 0, 2}), Monomial({3}), Monomial({2, 0, 2})
    };

    EXPECT_TRUE(GradedLexOrder{}(monomials[3], monomials[0]));
    EXPECT_TRUE(GradedLexOrder{}(monomials[0], monomials[2]));
    EXPECT_TRUE(GradedLexOrder{}(monomials[2], monomials[1]));

    EXPECT_TRUE(GradedLexOrder{}(monomials[3], monomials[2]));
    EXPECT_TRUE(GradedLexOrder{}(monomials[3], monomials[1]));
    EXPECT_TRUE(GradedLexOrder{}(monomials[0], monomials[1]));

    EXPECT_FALSE(GradedLexOrder{}(monomials[0], monomials[3]));
    EXPECT_FALSE(GradedLexOrder{}(monomials[2], monomials[0]));
    EXPECT_FALSE(GradedLexOrder{}(monomials[1], monomials[2]));

    EXPECT_FALSE(GradedLexOrder{}(monomials[2], monomials[3]));
    EXPECT_FALSE(GradedLexOrder{}(monomials[1], monomials[3]));
    EXPECT_FALSE(GradedLexOrder{}(monomials[1], monomials[0]));
}

TEST(MonomialOrderTest, GradedReverseLexOrderFull) {
    std::vector<Monomial> monomials = {
        Monomial({1, 2, 1}), Monomial({0, 0, 2}), Monomial({3}), Monomial({2, 0, 2})
    };

    EXPECT_TRUE(GradedReverseLexOrder{}(monomials[0], monomials[3]));
    EXPECT_TRUE(GradedReverseLexOrder{}(monomials[3], monomials[2]));
    EXPECT_TRUE(GradedReverseLexOrder{}(monomials[2], monomials[1]));

    EXPECT_TRUE(GradedReverseLexOrder{}(monomials[0], monomials[2]));
    EXPECT_TRUE(GradedReverseLexOrder{}(monomials[0], monomials[1]));
    EXPECT_TRUE(GradedReverseLexOrder{}(monomials[3], monomials[1]));

    EXPECT_FALSE(GradedReverseLexOrder{}(monomials[3], monomials[0]));
    EXPECT_FALSE(GradedReverseLexOrder{}(monomials[2], monomials[3]));
    EXPECT_FALSE(GradedReverseLexOrder{}(monomials[1], monomials[2]));

    EXPECT_FALSE(GradedReverseLexOrder{}(monomials[2], monomials[0]));
    EXPECT_FALSE(GradedReverseLexOrder{}(monomials[1], monomials[0]));
    EXPECT_FALSE(GradedReverseLexOrder{}(monomials[1], monomials[3]));
}

TEST(MonomialOrderTest, TextbookPolynomialExampleHasExpectedOrdering) {
    const Monomial x_y_squared_z{1, 2, 1};
    const Monomial z_squared{0, 0, 2};
    const Monomial x_cubed{3};
    const Monomial x_squared_z_squared{2, 0, 2};
    const std::vector<Monomial> input{
        x_y_squared_z,
        z_squared,
        x_cubed,
        x_squared_z_squared,
    };

    auto lexicographic = input;
    std::sort(lexicographic.begin(), lexicographic.end(), LexicographicOrder{});
    EXPECT_EQ(
        lexicographic, (std::vector<Monomial>{
                           x_cubed,
                           x_squared_z_squared,
                           x_y_squared_z,
                           z_squared,
                       })
    );

    auto graded_lexicographic = input;
    std::sort(graded_lexicographic.begin(), graded_lexicographic.end(), GradedLexOrder{});
    EXPECT_EQ(
        graded_lexicographic, (std::vector<Monomial>{
                                  x_squared_z_squared,
                                  x_y_squared_z,
                                  x_cubed,
                                  z_squared,
                              })
    );

    auto graded_reverse_lexicographic = input;
    std::sort(
        graded_reverse_lexicographic.begin(), graded_reverse_lexicographic.end(),
        GradedReverseLexOrder{}
    );
    EXPECT_EQ(
        graded_reverse_lexicographic, (std::vector<Monomial>{
                                          x_y_squared_z,
                                          x_squared_z_squared,
                                          x_cubed,
                                          z_squared,
                                      })
    );
}

TEST(MonomialOrderTest, EqualCanonicalMonomialsAreEquivalentInEveryOrder) {
    const Monomial lhs{2, 0, 3, 0, 0};
    const Monomial rhs{2, 0, 3};

    EXPECT_FALSE(LexicographicOrder{}(lhs, rhs));
    EXPECT_FALSE(LexicographicOrder{}(rhs, lhs));
    EXPECT_FALSE(GradedLexOrder{}(lhs, rhs));
    EXPECT_FALSE(GradedLexOrder{}(rhs, lhs));
    EXPECT_FALSE(GradedReverseLexOrder{}(lhs, rhs));
    EXPECT_FALSE(GradedReverseLexOrder{}(rhs, lhs));
}

TEST(MonomialOrderTest, GradedOrdersPrioritizeTotalDegree) {
    const Monomial lower_degree{100};
    const Monomial higher_degree{0, 0, 101};

    EXPECT_TRUE(GradedLexOrder{}(higher_degree, lower_degree));
    EXPECT_FALSE(GradedLexOrder{}(lower_degree, higher_degree));
    EXPECT_TRUE(GradedReverseLexOrder{}(higher_degree, lower_degree));
    EXPECT_FALSE(GradedReverseLexOrder{}(lower_degree, higher_degree));
}

TEST(MonomialOrderDetailTest, TotalDegreeOrderTreatsEqualDegreesAsEquivalent) {
    const Monomial first{4};
    const Monomial same_degree{1, 3};
    const Monomial higher_degree{0, 0, 5};

    EXPECT_FALSE(detail::TotalDegreeOrder{}(first, same_degree));
    EXPECT_FALSE(detail::TotalDegreeOrder{}(same_degree, first));
    EXPECT_TRUE(detail::TotalDegreeOrder{}(higher_degree, first));
    EXPECT_FALSE(detail::TotalDegreeOrder{}(first, higher_degree));
}

TEST(MonomialOrderDetailTest, ReverseLexOrderUsesTheLastDifferingVariable) {
    const Monomial x_squared{2};
    const Monomial x_y{1, 1};
    const Monomial equal_with_trailing_zeros{2, 0, 0};

    EXPECT_TRUE(detail::ReverseLexOrder{}(x_squared, x_y));
    EXPECT_FALSE(detail::ReverseLexOrder{}(x_y, x_squared));
    EXPECT_FALSE(detail::ReverseLexOrder{}(x_squared, equal_with_trailing_zeros));
    EXPECT_FALSE(detail::ReverseLexOrder{}(equal_with_trailing_zeros, x_squared));
}

TEST(MonomialOrderDetailTest, CompositeOrderUsesSecondaryOrderOnlyOnPrimaryTies) {
    using DegreeThenReverseLex =
        detail::CompositeOrder<detail::TotalDegreeOrder, detail::ReverseLexOrder>;

    const Monomial lower_degree{100};
    const Monomial higher_degree{0, 0, 101};
    const Monomial x_squared{2};
    const Monomial x_y{1, 1};

    EXPECT_TRUE(DegreeThenReverseLex{}(higher_degree, lower_degree));
    EXPECT_TRUE(DegreeThenReverseLex{}(x_squared, x_y));
    EXPECT_FALSE(DegreeThenReverseLex{}(x_y, x_squared));
}

TEST(MonomialOrderTest, OrdersProduceExpectedSortedSequences) {
    const std::vector<Monomial> input{
        Monomial{}, Monomial{0, 1}, Monomial{1}, Monomial{0, 2}, Monomial{1, 1}, Monomial{2},
    };

    auto lexicographic = input;
    std::sort(lexicographic.begin(), lexicographic.end(), LexicographicOrder{});
    EXPECT_EQ(
        lexicographic, (std::vector<Monomial>{
                           Monomial{2},
                           Monomial{1, 1},
                           Monomial{1},
                           Monomial{0, 2},
                           Monomial{0, 1},
                           Monomial{},
                       })
    );

    auto graded_lexicographic = input;
    std::sort(graded_lexicographic.begin(), graded_lexicographic.end(), GradedLexOrder{});
    EXPECT_EQ(
        graded_lexicographic, (std::vector<Monomial>{
                                  Monomial{2},
                                  Monomial{1, 1},
                                  Monomial{0, 2},
                                  Monomial{1},
                                  Monomial{0, 1},
                                  Monomial{},
                              })
    );

    auto graded_reverse_lexicographic = input;
    std::sort(
        graded_reverse_lexicographic.begin(), graded_reverse_lexicographic.end(),
        GradedReverseLexOrder{}
    );
    EXPECT_EQ(graded_reverse_lexicographic, graded_lexicographic);
}

namespace {

template <typename Order>
void ExpectStrictWeakOrdering(const std::vector<Monomial>& monomials) {
    const Order order;

    for (const auto& lhs : monomials) {
        EXPECT_FALSE(order(lhs, lhs));
        for (const auto& rhs : monomials) {
            if (order(lhs, rhs)) {
                EXPECT_FALSE(order(rhs, lhs));
            }
            for (const auto& third : monomials) {
                if (order(lhs, rhs) && order(rhs, third)) {
                    EXPECT_TRUE(order(lhs, third));
                }
            }
        }
    }
}

}  // namespace

TEST(MonomialOrderTest, EveryOrderSatisfiesStrictWeakOrderingLaws) {
    const std::vector<Monomial> monomials{
        Monomial{},     Monomial{1},    Monomial{0, 1},    Monomial{2},
        Monomial{1, 1}, Monomial{0, 2}, Monomial{1, 0, 1},
    };

    ExpectStrictWeakOrdering<LexicographicOrder>(monomials);
    ExpectStrictWeakOrdering<GradedLexOrder>(monomials);
    ExpectStrictWeakOrdering<GradedReverseLexOrder>(monomials);
    ExpectStrictWeakOrdering<detail::TotalDegreeOrder>(monomials);
    ExpectStrictWeakOrdering<detail::ReverseLexOrder>(monomials);
}
