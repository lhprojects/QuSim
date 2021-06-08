cmake_minimum_required (VERSION 3.8)

file(GLOB SRC QuSolverBenchmark.cpp)
message(STATUS "files: ${SRC}")


add_executable(QuSolver1DBenchmark ${SRC})

target_link_libraries(QuSolver1DBenchmark QuSim)
