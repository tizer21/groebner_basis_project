cmake_minimum_required(VERSION 3.20)

project(
    groebner_basis
    VERSION 0.1
    DESCRIPTION "Header-only Gröbner-basis library"
    LANGUAGES CXX
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(groebner_basis INTERFACE)

add_library(
    groebner_basis::groebner_basis
    ALIAS groebner_basis
)

target_compile_features(
    groebner_basis
    INTERFACE
        cxx_std_20
)

target_include_directories(
    groebner_basis
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

include(CTest)
include(FetchContent)

if(BUILD_TESTING)
    FetchContent_Declare(
        googletest
        URL
            https://github.com/google/googletest/archive/refs/tags/v1.18.0.zip
    )

    if(WIN32)
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    endif()

    FetchContent_MakeAvailable(googletest)
    add_subdirectory(tests)
endif()
