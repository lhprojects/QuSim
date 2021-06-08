cmake_minimum_required (VERSION 3.8)

file(GLOB QuScattering3DPrecisionSources QuScattering3DPrecision.cpp)
message(STATUS "files: ${QuScattering3DPrecisionSources}")


add_executable(QuScattering3DPrecision ${QuScattering3DPrecisionSources})

target_link_libraries(QuScattering3DPrecision QuSim)
