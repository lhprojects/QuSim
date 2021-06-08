cmake_minimum_required (VERSION 3.8)

file(GLOB QuSimSources QuSim.cpp QuSim.rc)

message(STATUS "files: ${QuSimSources}")


add_executable(QuEvol1D WIN32 ${QuSimSources})

target_link_libraries(QuEvol1D QuSim)
