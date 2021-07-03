cmake_minimum_required (VERSION 3.18)
project ("QuSim")

option(USE_CUDA "use cuda" OFF)
if(USE_CUDA)
    add_compile_definitions(QUSIM_USE_CUDA)
endif()

message("Compiler: ${CMAKE_CXX_COMPILER_ID}")

set (CMAKE_CXX_STANDARD 17)

if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

    # this can't deal with cuda host compiler correctly
    
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")


    ADD_DEFINITIONS(-DUNICODE)
    ADD_DEFINITIONS(-D_UNICODE)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-Wfatal-errors -g -march=native)
else()
    add_compile_options(-Wfatal-errors -g -march=native -stdlib=libc++)
    add_link_options(-stdlib=libc++ -lc++abi -lc++ -lc++abi -lm -lc -lgcc_s -lgcc)
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../third_party) # for eigen

set(HEADER_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/*.h
)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory("nana")

add_subdirectory("QuSimCore")
add_subdirectory("QuSim")
add_subdirectory("QuSimTest")
add_subdirectory("QuBenchmarkLoop")

add_subdirectory("QuCalBenchmark")

add_subdirectory("QuEvol1D")
add_subdirectory("QuEvol1DBenchmark")
add_subdirectory("QuEvol1DPrecision")

add_subdirectory("QuEvol2D")
add_subdirectory("QuEvol2DBenchmark")

add_subdirectory("QuScattering1D")
add_subdirectory("QuScattering1DPrecision")

add_subdirectory("QuScattering2D")
add_subdirectory("QuScattering2DPrecision")

add_subdirectory("QuScattering3DPrecision")

add_subdirectory("QuSolver1D")
add_subdirectory("QuSolver1DBenchmark")
add_subdirectory("QuSolver1DPrecision")


