cmake_minimum_required (VERSION 3.8)

file(GLOB QuSim2DBenchmarkSources QuSim2DBenchmark.cpp)
message(STATUS "files: ${QuSim2DBenchmarkSources}")


add_executable(QuEvol2DBenchmark ${QuSim2DBenchmarkSources})

target_link_libraries(QuEvol2DBenchmark QuSim)
