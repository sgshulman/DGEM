cmake_minimum_required(VERSION 3.5)
set(CMAKE_GENERATOR "Unix Makefiles")

project(DGEM)

set(SOURCES
        directions.cpp
        directions.hpp
        Direction3d.hpp
        grid.cpp
        grid.hpp
        observers.hpp
        model.cpp
        model.hpp
        photons.cpp
        photons.hpp
        points.cpp
        Vector3d.hpp)

add_executable(dgem ${SOURCES})

set_target_properties(dgem PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wextra -Wpedantic")


