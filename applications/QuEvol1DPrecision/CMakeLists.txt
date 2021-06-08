cmake_minimum_required (VERSION 3.8)

file(GLOB QuSimPrecisionSources QuSimPrecision.cpp)
message(STATUS "files: ${QuSimPrecisionSources}")


add_executable(QuEvol1DPrecision ${QuSimPrecisionSources})

target_link_libraries(QuEvol1DPrecision QuSim)
