cmake_minimum_required(VERSION 3.5)
set(CMAKE_GENERATOR "Unix Makefiles")

project(DGEM)

set(SOURCES
        AzimuthalHump.cpp
        AzimuthalHump.hpp
        CartesianGrid.cpp
        CartesianGrid.hpp
        DebugUtils.cpp
        DebugUtils.hpp
        Directions.cpp
        Directions.hpp
        Direction3d.cpp
        Direction3d.hpp
        FlaredDisk.cpp
        FlaredDisk.hpp
        FractalCloud.cpp
        FractalCloud.hpp
        Halton.cpp
        Halton.hpp
        IDiskHump.hpp
        IDust.hpp
        IGrid.hpp
        IMatter.hpp
        IRandomGenerator.hpp
        KurosawaWind.cpp
        KurosawaWind.hpp
        LEcuyer.cpp
        LEcuyer.hpp
        MathUtils.hpp
        MatterArray.cpp
        MatterArray.hpp
        MatterTranslation.cpp
        MatterTranslation.hpp
        MieDust.cpp
        MieDust.hpp
        model.cpp
        model.hpp
        Observer.cpp
        Observer.hpp
        Photon.cpp
        Photon.hpp
        Predefines.hpp
        RoundHump.cpp
        RoundHump.hpp
        Sources.cpp
        Sources.hpp
        SafierWind.cpp
        SafierWind.hpp
        Sobol.cpp
        Sobol.hpp
        SphereEnvelope.cpp
        SphereEnvelope.hpp
        StdRandomGenerator.cpp
        StdRandomGenerator.hpp
        TetrahedralGrid.cpp
        TetrahedralGrid.hpp
        Units.hpp
        Vector2d.hpp
        Vector3d.hpp
        WhiteDust.cpp
        WhiteDust.hpp
        third-party/nlohmann/json.hpp)

if (DEFINED CLANG_TIDY_VERSION)
    message("CLANG_TIDY ENABLED: ${CLANG_TIDY_VERSION}")
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_VERSION}")
endif()

add_executable(dgem points.cpp ${SOURCES})

set_target_properties(dgem PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO)

if(DEFINED ENABLE_UNIT_TESTS AND ${ENABLE_UNIT_TESTS})
    message("UNIT TESTS ENABLED")

    set(UNITTEST_SOURCES
        tests/unittests.cpp
        tests/AzimuthalHumpTests.cpp
        tests/CartesianGridTests.cpp
        tests/Direction3dTests.cpp
        tests/HaltonTests.cpp
        tests/LEcuyerTests.cpp
        tests/MathUtilsTests.cpp
        tests/MatterArrayTests.cpp
        tests/MatterTranslationTests.cpp
        tests/MieDustTests.cpp
        tests/ModelTests.inl
        tests/ObserverTests.cpp
        tests/RoundHumpTests.cpp
        tests/SobolTests.cpp
        tests/SourcesTests.cpp
        tests/SphereEnvelopeTests.cpp
        tests/StdRandomGeneratorTests.cpp
        tests/TestUtils.hpp
        tests/TetrahedralGridTests.cpp
        tests/Vector3dTests.cpp
        tests/WhiteDustTests.cpp
        third-party/catch2/catch.hpp)

    # MacOs 11 Catalina has problems with cmath
    if (NOT (APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
        set(UNITTEST_SOURCES ${UNITTEST_SOURCES} tests/DirectionsTests.cpp)
    endif()

    add_executable(unittests ${UNITTEST_SOURCES} ${SOURCES})

    set_target_properties(unittests PROPERTIES
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
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -static-libstdc++ -Wa,-mbig-obj")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -static-libstdc++ -Wa,-mbig-obj")
endif()
