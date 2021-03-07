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
        Direction3d.hpp
        FlaredDisk.cpp
        FlaredDisk.hpp
        FractalCloud.cpp
        FractalCloud.hpp
        IDiskHump.hpp
        IDust.hpp
        IGrid.hpp
        IMatter.hpp
        KurosawaWind.cpp
        KurosawaWind.hpp
        observers.hpp
        MathUtils.hpp
        MatterArray.cpp
        MatterArray.hpp
        MatterTranslation.cpp
        MatterTranslation.hpp
        MieDust.cpp
        MieDust.hpp
        model.cpp
        model.hpp
        Photon.cpp
        Photon.hpp
        Predefines.hpp
        Random.hpp
        RoundHump.cpp
        RoundHump.hpp
        Sources.cpp
        Sources.hpp
        SafierWind.cpp
        SafierWind.hpp
        SphereEnvelope.cpp
        SphereEnvelope.hpp
        TetrahedralGrid.cpp
        TetrahedralGrid.hpp
        Units.hpp
        Vector3d.hpp
        WhiteDust.cpp
        WhiteDust.hpp
        third-party/nlohmann/json.hpp)

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
        tests/MathUtilsTests.cpp
        tests/MatterArrayTests.cpp
        tests/MatterTranslationTests.cpp
        tests/MieDustTests.cpp
        tests/ModelTests.inl
        tests/RoundHumpTests.cpp
        tests/SphereEnvelopeTests.cpp
        tests/TestUtils.hpp
        tests/Vector3dTests.cpp
        tests/WhiteDustTests.cpp
        third-party/catch2/catch.hpp)

    add_executable(unittests ${UNITTEST_SOURCES} ${SOURCES})

    set_target_properties(unittests PROPERTIES
            CXX_STANDARD 14
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
