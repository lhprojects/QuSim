cmake_minimum_required (VERSION 3.8)

file(GLOB SRC QuSolverPrecision.cpp)
message(STATUS "files: ${SRC}")


add_executable(QuSolver1DPrecision ${SRC})

target_link_libraries(QuSolver1DPrecision QuSim)
