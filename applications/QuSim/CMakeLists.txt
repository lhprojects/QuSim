cmake_minimum_required (VERSION 3.8)

file(GLOB SRC ../../src/QuSim.cpp)
add_compile_definitions(QUEXPORT)


message(STATUS "QuSim --         files: ${SRC}")
add_library(QuSim SHARED ${SRC})
target_link_libraries(QuSim QuSimCore)

message(STATUS "QuSim -- release flags: " ${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS "QuSim --   debug falgs: " ${CMAKE_CXX_FLAGS_DEBUG})
