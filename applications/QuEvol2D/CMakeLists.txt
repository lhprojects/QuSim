cmake_minimum_required (VERSION 3.8)

file(GLOB QuSim2DSources QuSim2D.cpp QuSim2D.rc)

message(STATUS "files: ${QuSim2DSources}")


add_executable(QuEvol2D WIN32 ${QuSim2DSources})

target_link_libraries(QuEvol2D QuSim)
