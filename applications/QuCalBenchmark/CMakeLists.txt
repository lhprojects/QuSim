cmake_minimum_required (VERSION 3.8)

file(GLOB SRC *.cpp)

message(STATUS "files: ${SRC}")


add_executable(QuCalBenchmark ${SRC})

target_link_libraries(QuCalBenchmark QuSim)
