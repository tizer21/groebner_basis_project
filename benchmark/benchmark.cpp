#include <benchmark/benchmark.h>
#include <stddef.h>

#include "generators.h"
#include "groebner_basis/algorithms.h"
#include "groebner_basis/fields/prime_field.h"
#include "groebner_basis/monomial_orders.h"

namespace {

namespace algorithms = groebner_basis::algorithms;
namespace benchmarking = groebner_basis::benchmarking;

using Coefficient = groebner_basis::fields::PrimeField<32003>;
using Improvement = algorithms::detail::Improvement;

using LexicographicOrder = groebner_basis::LexicographicOrder;
using GradedLexOrder = groebner_basis::GradedLexOrder;
using GradedReverseLexOrder = groebner_basis::GradedReverseLexOrder;

struct BasisSizes {
    size_t groebner;
    size_t minimal;
    size_t reduced;
};

struct ImprovementCase {
    std::string_view name;
    Improvement improvement;
};

constexpr std::array<ImprovementCase, 3> kImprovementCases{{
    {"noCriteria", Improvement::kNone},
    {"productCriterion", Improvement::kProductCriterion},
    {"productAndChainCriteria", Improvement::kProductAndChainCriteria},
}};

template <typename PolynomialBasisType>
void BenchmarkSystem(
    benchmark::State& state, const PolynomialBasisType& generators, Improvement improvement
) {
    BasisSizes sizes;
    for (auto _ : state) {
        state.PauseTiming();
        auto input = generators;
        state.ResumeTiming();

        auto groebner_basis =
            algorithms::detail::RunBuchbergerAlgorithm(std::move(input), improvement);
        sizes.groebner = groebner_basis.size();

        auto minimal_basis = algorithms::detail::MinimizeGroebnerBasis(std::move(groebner_basis));
        sizes.minimal = minimal_basis.size();

        auto reduced_basis = algorithms::detail::InterreduceGroebnerBasis(std::move(minimal_basis));
        sizes.reduced = reduced_basis.size();

        benchmark::DoNotOptimize(reduced_basis);
        benchmark::ClobberMemory();
    }

    state.counters["groebner_basis_size"] = static_cast<double>(sizes.groebner);
    state.counters["minimal_basis_size"] = static_cast<double>(sizes.minimal);
    state.counters["reduced_basis_size"] = static_cast<double>(sizes.reduced);
}

template <typename MonomialOrder>
void BenchmarkCyclic(benchmark::State& state, Improvement improvement) {
    const auto variable_count = static_cast<std::size_t>(state.range(0));
    const auto generators =
        benchmarking::GenerateCyclicSystem<Coefficient, MonomialOrder>(variable_count);

    BenchmarkSystem(state, generators, improvement);
}

template <typename MonomialOrder>
void BenchmarkKatsura(benchmark::State& state, Improvement improvement) {
    const auto order = static_cast<std::size_t>(state.range(0));
    const auto generators = benchmarking::GenerateKatsuraSystem<Coefficient, MonomialOrder>(order);

    BenchmarkSystem(state, generators, improvement);
}

std::string MakeBenchmarkName(
    std::string_view system, std::string_view monomial_order, std::string_view improvement
) {
    std::string name{"ReducedBasis/"};

    name.append(system);
    name.push_back('/');
    name.append(monomial_order);
    name.push_back('/');
    name.append(improvement);

    return name;
}

template <typename MonomialOrder>
void RegisterBenchmarksForOrder(std::string_view order_name) {
    for (const auto& benchmark_case : kImprovementCases) {
        auto* cyclic_benchmark = benchmark::RegisterBenchmark(
            MakeBenchmarkName("Cyclic", order_name, benchmark_case.name),
            [improvement = benchmark_case.improvement](benchmark::State& state) {
                BenchmarkCyclic<MonomialOrder>(state, improvement);
            }
        );

        cyclic_benchmark->DenseRange(3, 4)->ArgName("variables")->Unit(benchmark::kSecond);

        auto* katsura_benchmark = benchmark::RegisterBenchmark(
            MakeBenchmarkName("Katsura", order_name, benchmark_case.name),
            [improvement = benchmark_case.improvement](benchmark::State& state) {
                BenchmarkKatsura<MonomialOrder>(state, improvement);
            }
        );
        katsura_benchmark->DenseRange(1, 3)->ArgName("order")->Unit(benchmark::kSecond);
    }
}

const bool kBenchmarksRegistered = [] {
    RegisterBenchmarksForOrder<groebner_basis::LexicographicOrder>("Lexicographic");
    RegisterBenchmarksForOrder<groebner_basis::GradedLexOrder>("GradedLex");
    RegisterBenchmarksForOrder<groebner_basis::GradedReverseLexOrder>("GradedReverseLex");

    return true;
}();

}  // namespace
