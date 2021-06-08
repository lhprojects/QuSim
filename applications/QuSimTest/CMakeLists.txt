cmake_minimum_required (VERSION 3.8)

file(GLOB QuSimTestSources Test.cpp)

set(QuSimTestSources
${QuSimTestSources}
${CMAKE_CURRENT_SOURCE_DIR}/../../src/TestCal.cpp)

message(STATUS "files: ${QuSimTestSources}")

add_executable(QuSimTest ${QuSimTestSources})
target_link_libraries(QuSimTest QuSimCore)
