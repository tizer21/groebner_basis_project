# Программный проект «Базисы Грёбнера»

Библиотека на C++20 для работы с многочленами от нескольких переменных,
идеалами над кольцом многочленов и базисами Грёбнера.

## Возможности

- арифметика рациональных чисел и конечных простых полей;
- представление мономов, термов и многочленов от нескольких переменных;
- лексикографический (`Lexicographic`), однородный лексикографический (`GradedLex`) и
  однородный обратный лексикографический (`GradedRevLex`) порядки;
- деление многочленов, редукция и вычисление S-многочленов;
- алгоритм Бухбергера с возможностью выбора метода оптимизации (критерии произведения и цепочки);
- вычисление обычного, минимального и редуцированного базисов Грёбнера;
- проверка принадлежности многочлена идеалу и вычисление нормальной формы;
- сложение и умножение идеалов над кольцом многочленов;
- юнит-тесты и набор бенчмарков.

## Требования

- компилятор с поддержкой C++20;
- CMake 3.20 или новее;
- доступ к интернету при первой сборке тестов или бенчмарков - CMake загружает
  GoogleTest или Google Benchmark через `FetchContent`.

Сама библиотека является заголовочной и не имеет внешних зависимостей.

## Пример использования

Вектор степеней в `Monomial` задаётся в порядке
`{степень x_0, степень x_1, ...}`. Например, `Monomial{2, 1}` соответствует
моному `x_0^2 * x_1`.

```cpp
#include <iostream>

#include "groebner_basis/algorithms.h"
#include "groebner_basis/fields/rational.h"
#include "groebner_basis/monomial_orders.h"

int main() {
    using namespace groebner_basis;

    using Coefficient = fields::Rational;
    using TermType = Term<Coefficient>;
    using PolynomialType = Polynomial<Coefficient, LexicographicOrder>;
    using Basis = algorithms::PolynomialBasis<Coefficient, LexicographicOrder>;

    // f_1 = x_0^2 - x_1
    const PolynomialType first{
        TermType{Coefficient{1}, Monomial{2}},
        TermType{Coefficient{-1}, Monomial{0, 1}},
    };

    // f_2 = x_0*x_1 - 1
    const PolynomialType second{
        TermType{Coefficient{1}, Monomial{1, 1}},
        TermType{Coefficient{-1}},
    };

    const auto basis =
        algorithms::ComputeReducedGroebnerBasis(Basis{first, second});

    for (const auto& polynomial : basis) {
        std::cout << polynomial << '\n';
    }
}
```

Для подключения библиотеки к другому CMake-проекту:

```cmake
add_subdirectory(path/to/groebner_basis_project)
target_link_libraries(your_target PRIVATE groebner_basis::groebner_basis)
```

## Сборка и запуск тестов

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Тесты написаны с использованием GoogleTest.

## Бенчмарки

Бенчмарки сравнивают варианты алгоритма Бухбергера на системах Cyclic и
Katsura для поддерживаемых порядков мономов.

```bash
cmake -S . -B build-benchmark \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF \
  -DGROEBNER_BASIS_BUILD_BENCHMARKS=ON

cmake --build build-benchmark --target groebner_basis_benchmarks -j
./build-benchmark/benchmark/groebner_basis_benchmarks
```

Список доступных бенчмарков:

```bash
./build-benchmark/benchmark/groebner_basis_benchmarks --benchmark_list_tests
```

Бенчмарки предназначены для измерения производительности, а не для проверки
корректности.

## Структура проекта

```text
include/groebner_basis/  публичные заголовочные файлы
tests/                   модульные тесты
benchmark/               генераторы и бенчмарки
```
