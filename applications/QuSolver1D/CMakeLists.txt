cmake_minimum_required (VERSION 3.8)

file(GLOB QuSolver1DSource QuSolver.cpp QuSolver.rc)

message(STATUS "files: ${QuSolver1DSource}")


add_executable(QuSolver1D WIN32 ${QuSolver1DSource})

target_link_libraries(QuSolver1D QuSim)
