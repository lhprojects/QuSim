cmake_minimum_required (VERSION 3.8)

file(GLOB SRC QuBenchmarkLoop.cpp)
message(STATUS "files: ${SRC}")


add_executable(QuBenchmarkLoop ${SRC})

target_link_libraries(QuBenchmarkLoop QuSimCore)

