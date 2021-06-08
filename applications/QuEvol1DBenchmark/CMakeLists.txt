cmake_minimum_required (VERSION 3.8)

file(GLOB QuSimBenchmarkSources QuSimBenchmark.cpp)
message(STATUS "files: ${QuSimBenchmarkSources}")


add_executable(QuEvol1DBenchmark ${QuSimBenchmarkSources})

target_link_libraries(QuEvol1DBenchmark QuSim)
