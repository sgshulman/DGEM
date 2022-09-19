cmake_minimum_required(VERSION 3.5)
set(CMAKE_GENERATOR "Unix Makefiles")

project(UTILS)

set(IMDIFF_SOURCES
        Image.hpp
        imdiff.cpp)

add_executable(imdiff ${IMDIFF_SOURCES})

set_target_properties(imdiff PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO)

set(IMSUM_SOURCES
        Image.hpp
        imsum.cpp)

add_executable(imsum ${IMSUM_SOURCES})

set_target_properties(imsum PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO)

if(DEFINED ENABLE_UNIT_TESTS AND ${ENABLE_UNIT_TESTS})
    message("UNIT TESTS ENABLED")

    set(UNITTEST_SOURCES
        Image.hpp
        imdiff.cpp
        imsum.cpp
        tests/DifferTests.inl
        tests/SumTests.inl
        tests/utilstests.cpp
        ../source/third-party/catch2/catch.hpp)

    add_executable(utiltests ${UNITTEST_SOURCES})

    set_target_properties(utiltests PROPERTIES
            CXX_STANDARD 17
            CXX_STANDARD_REQUIRED YES
            CXX_EXTENSIONS NO
            COMPILE_DEFINITIONS "ENABLE_UNIT_TESTS")
endif()

if ( MSVC )
    set(CMAKE_CXX_FLAGS_DEBUG "/W4 /WX")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /W4 /WX")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -Wpedantic")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wextra -Wpedantic")
endif()

# MinGW has problems with the linker
if ( MINGW )
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -static-libstdc++")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -static-libstdc++")
endif()
