cmake_minimum_required(VERSION 3.22)
project(Assignment1_2024)

set(CMAKE_CXX_STANDARD 23)

set(HEADERS
        headers/matrixhelpers.h
        headers/squarematrix.h
        headers/vectorhelpers.h
        headers/shift_inverse_power_iteration.h
        headers/eig_finder_helpers.h
        headers/power_iteration.h
        headers/inverse_power_iteration.h
)

set(SOURCES
        src/matrixhelpers.cpp
        src/squarematrix.cpp
        src/vectorhelpers.cpp
        src/power_iteration.cpp
        src/inverse_power_iteration.cpp
        src/shift_inverse_power_iteration.cpp
)

add_executable(Assignment1_2024 main.cpp ${HEADERS} ${SOURCES})