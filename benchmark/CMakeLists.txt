add_executable(
    groebner_basis_benchmarks
    benchmark.cpp
)

target_link_libraries(
    groebner_basis_benchmarks
    PRIVATE
        groebner_basis::groebner_basis
        benchmark::benchmark
        benchmark::benchmark_main
)

target_include_directories(
    groebner_basis_benchmarks
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(
    groebner_basis_benchmarks
    PRIVATE
        cxx_std_20
)
