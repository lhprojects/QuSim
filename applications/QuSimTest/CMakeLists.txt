cmake_minimum_required (VERSION 3.8)

file(GLOB SRC Test.cpp TestFFT.cpp)

set(SRC 
${SRC}
${CMAKE_CURRENT_SOURCE_DIR}/../../src/TestCal.cpp)

message(STATUS "files: ${SRC}")

add_executable(QuSimTest ${SRC})
target_link_libraries(QuSimTest QuSimCore)
