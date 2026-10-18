add_executable(
    groebner_tests
    test_monomial.cpp
    test_prime_field.cpp
    test_rational.cpp
)

target_link_libraries(
    groebner_tests
    PRIVATE
        groebner_basis::groebner_basis
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(groebner_tests)