cmake_minimum_required (VERSION 3.17)

SET(SRC
Cal.cpp
ComplexPotentialIVPSolver1DImpl.cpp
EigenMethod.cpp
EvolverImpl.cpp
FourierTransform.cpp
GaussLegendreMethod1D.cpp
GaussLegendreMethod2D.cpp
IVPSolverImpl.cpp
MatrixSolver.cpp
OptionsImpl.cpp
Perturbation1DImpl.cpp
Perturbation2DImpl.cpp
Perturbation3DImpl.cpp
PerturbationOptions.cpp
ScatteringProblemSolverInverseMatrix1D.cpp
ScatteringProblemSolverInverseMatrix2D.cpp
ScatteringProblemSolverInverseMatrix3D.cpp
ScatteringSolverImpl.cpp
SolverImpl.cpp
SplittingMethod1D.cpp
SplittingMethod1DCUDA.cpp
SplittingMethod2D.cpp
SplittingMethod2DCUDA.cpp
)

option(USE_CUDA "use cuda" OFF)
if(USE_CUDA)
    find_package(CUDAToolkit)
    find_package(CUDA REQUIRED)

    list(APPEND SRC CudaUtility.cu)
    #enable_language(CUDA)
    add_compile_definitions(QUSIM_USE_CUDA)
    add_compile_definitions(USE_CUDA)


else()
    message("DON'T USE CUDA")
endif()


list(TRANSFORM SRC PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/../../src/)
message(STATUS "files: ${SRC}")


if(USE_CUDA)
cuda_add_library(QuSimCore ${SRC})
cuda_add_cufft_to_target(QuSimCore)
cuda_add_cublas_to_target(QuSimCore)

#target_link_libraries(QuSimCore PRIVATE CUDA::cudart)
#target_link_libraries(QuSimCore PRIVATE CUDA::cublas)
#target_link_libraries(QuSimCore PRIVATE CUDA::cufft)

else()
add_library(QuSimCore ${SRC})
endif()

