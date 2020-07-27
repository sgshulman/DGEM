cmake_minimum_required(VERSION 3.5)
set(CMAKE_GENERATOR "Unix Makefiles")

project(DGEM)

set(SOURCES
        Directions.cpp
        Directions.hpp
        Direction3d.hpp
        Dust.cpp
        Dust.hpp
        FlaredDisk.cpp
        FlaredDisk.hpp
        grid.cpp
        grid.hpp
        IMatter.hpp
        observers.hpp
        MathUtils.hpp
        MatterArray.cpp
        MatterArray.hpp
        model.cpp
        model.hpp
        Photon.cpp
        Photon.hpp
        Random.hpp
        Sources.cpp
        Sources.hpp
        SafierWind.cpp
        SafierWind.hpp
        SphereEnvelope.cpp
        SphereEnvelope.hpp
        Predefines.hpp
        Vector3d.hpp
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
        tests/Direction3dTests.cpp
        tests/MathUtilsTests.cpp
        tests/MatterArrayTests.cpp
        tests/ModelTests.inl
        tests/SphereEnvelopeTests.cpp
        tests/Vector3dTests.cpp
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
