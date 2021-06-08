cmake_minimum_required (VERSION 3.8)

file(GLOB QuSimSources ../../src/QuSim.cpp)
add_compile_definitions(QUEXPORT)

message(STATUS "files: ${QuSimSources}")
add_library(QuSim SHARED ${QuSimSources})
target_link_libraries(QuSim QuSimCore)
