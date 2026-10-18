add_executable(
    groebner_tests
    test_algorithms.cpp
    test_monomial_orders.cpp
    test_monomial.cpp
    test_polynomial_ideal.cpp
    test_polynomial.cpp
    test_prime_field.cpp
    test_rational.cpp
    test_term.cpp
)

target_link_libraries(
    groebner_tests
    PRIVATE
        groebner_basis::groebner_basis
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(groebner_tests)