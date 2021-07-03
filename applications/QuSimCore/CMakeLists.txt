cmake_minimum_required (VERSION 3.17)

option(USE_OPENMP "use OpenMP" OFF)


SET(SRC
Cal.cpp
ComplexPotentialIVPSolver1DImpl.cpp
Device.cpp
EigenMethod1D.cpp
EvolverImpl.cpp
EvolverImpl1D.cpp
EvolverImpl2D.cpp
FourierTransform.cpp
GaussLegendreMethod.cpp
GaussLegendreMethod1D.cpp
GaussLegendreMethod2D.cpp
IVPSolverImpl.cpp
MatrixSolver.cpp
OptionsImpl.cpp
Perturbation1DImpl.cpp
Perturbation2DImpl.cpp
Perturbation3DImpl.cpp
PerturbationOptions.cpp
ScatteringProblemSolverInverseMatrix.cpp
ScatteringProblemSolverInverseMatrix1D.cpp
ScatteringProblemSolverInverseMatrix2D.cpp
ScatteringProblemSolverInverseMatrix3D.cpp
ScatteringSolverImpl.cpp
SolverImpl.cpp
SplittingUtils.cpp
SplittingMethod1D.cpp
SplittingMethod2D.cpp
)

file(GLOB HDR  ${CMAKE_CURRENT_SOURCE_DIR}/../../src/*.h)

option(USE_CUDA "use cuda" OFF)
if(USE_CUDA)
    find_package(CUDAToolkit)
    find_package(CUDA REQUIRED)

    set(CUDA_USE_STATIC_CUDA_RUNTIME ON)

    #enable_language(CUDA)

    list(APPEND SRC CudaDevice.cu)
    add_compile_definitions(QUSIM_USE_CUDA)
else()
    message("DON'T USE CUDA")
endif()


list(TRANSFORM SRC PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/../../src/)
message(STATUS "files: ${SRC}")
  
if(USE_CUDA)

    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-std=c++17)
    cuda_add_library(QuSimCore ${SRC} ${HDR})
    cuda_add_cufft_to_target(QuSimCore)
    cuda_add_cufft_to_target(QuSimCore)
    cuda_add_cublas_to_target(QuSimCore)

    #target_link_libraries(QuSimCore PRIVATE CUDA::cudart)
    #target_link_libraries(QuSimCore PRIVATE CUDA::cublas)
    #target_link_libraries(QuSimCore PRIVATE CUDA::cufft)

    message(STATUS ${CMAKE_MSVC_RUNTIME_LIBRARY})

    message(STATUS "QuSimCore -- CUDA_NVCC_FLAGS: " ${CUDA_NVCC_FLAGS})
    message(STATUS "QuSimCore -- CUDA_NVCC_FLAGS_DEBUG: " ${CUDA_NVCC_FLAGS_DEBUG})
    message(STATUS "QuSimCore -- CUDA_NVCC_FLAGS_RELEASE: " ${CUDA_NVCC_FLAGS_RELEASE})

    message(STATUS "QuSimCore -- CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
    message(STATUS "QuSimCore -- CMAKE_CXX_FLAGS_DEBUG: " ${CMAKE_CXX_FLAGS_DEBUG})
    message(STATUS "QuSimCore -- CMAKE_CXX_FLAGS_RELEASE: " ${CMAKE_CXX_FLAGS_RELEASE})

else()
    add_library(QuSimCore ${SRC} ${HDR})
endif()

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(QuSimCore OpenMP::OpenMP_CXX)
endif()

# Copy DLL
if(MSVC)
target_link_libraries(QuSimCore  ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/fftw-3.3.5-dll64/libfftw3-3.lib)

add_custom_target(CopyDll ALL)
add_custom_command(TARGET CopyDll PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/fftw-3.3.5-dll64/libfftw3-3.dll ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
)
endif()

