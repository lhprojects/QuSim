cmake_minimum_required (VERSION 3.8)

file(GLOB SRC QuScattering2D.cpp QuScattering2D.rc)

message(STATUS "files: ${SRC}")

add_executable(QuScattering2D WIN32 ${SRC})

target_link_libraries(QuScattering2D QuSim)
