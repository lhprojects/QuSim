cmake_minimum_required (VERSION 3.8)

file(GLOB QuScattering1DSources QuScattering1D.cpp)

message(STATUS "files: ${QuScattering1DSources}")


add_executable(QuScattering1D WIN32 ${QuScattering1DSources})

target_link_libraries(QuScattering1D QuSim)
target_link_libraries(QuScattering1D nana)
