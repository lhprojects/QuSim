cmake_minimum_required (VERSION 3.8)

file(GLOB QuScattering2DPrecisionSources QuScattering2DPrecision.cpp)
message(STATUS "files: ${QuScattering2DPrecisionSources}")


add_executable(QuScattering2DPrecision ${QuScattering2DPrecisionSources})

target_link_libraries(QuScattering2DPrecision QuSim)
