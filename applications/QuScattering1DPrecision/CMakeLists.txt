cmake_minimum_required (VERSION 3.8)

file(GLOB QuScattering1DPrecisionSources QuScattering1DPrecision.cpp)
message(STATUS "files: ${QuScattering1DPrecisionSources}")


add_executable(QuScattering1DPrecision ${QuScattering1DPrecisionSources})

target_link_libraries(QuScattering1DPrecision QuSim)
