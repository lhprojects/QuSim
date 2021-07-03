# QuSim
Quantum simulation and visualization program for wave function.

## Note
The core functions of this project is avoiding using the Windows SDK, except the GUI part.
Now the GUI part of this project is moving from native Windows API to the partable GUI library, ```nana```.
## Algorithm and Theory
 See the [FILE](theory/theory.pdf). Now the progress of document is far behind the codes.

## Quantum Evolution for 1 Dimension Wave Function
Solve the 1 D time evolution problem for given initial wave function and potential.
![Quantum Evolution for 1 Dimsional Wave Function](screenshot/QuSim1DExample.png)

## Quantum Evolution for 2 Dimension Wave Function
Solve the 2 D time evolution problem for given inital wave function and potential.
![Quantum Evolution for 2 Dimsional Wave Function](screenshot/QuSim2DExample.png)

## Solver for Quantum 1 Dimension Initial Value Problem
Solve the 1 D initial value problem for given intial condition, energy, and potential.
The solver for two-side boundary value problem using `shooting` method is also provoided. 
![Quantum Evolution for 1 Dimsional Wave Function](screenshot/QuSolverExample.png)

## Solver for Quantum 2 Dimension Scattering Problem
Solve the 2 D scattering problem for given initial wave function, energy, and potential.
![Quantum Scattering for 2 Dimsional Wave Function](screenshot/QuScattering2DExample.png)



## Q & A

Q: The difference between FFTW3 and KISS

A:  The FFTW3 is about as twice fast as the KISS.



## Hardware

GPU is supported by GPU_CUDA device. You can utilize GPU by building with CUDA and set options with `Options opts; opts.Cuda();`.
Multiple-threading is supported by CPU_PAR device. You can utilize multiple-core with `Options opts; opts.CpuPar();`.

Vectorized-instruction is supported by automatic vectorization support of the compiler. The `CpuParVec()` means we will try to use `std::execution::parallel_unsequenced_policy`. But it seems the same as the `std::execution::parallel_policy` . Similar, the `CpuVec()`  is the same as the `CpuSeq()`. However,  for `CpuParVec()` , `CpuVec()` , `CpuPar()` and `CpuSeq()`, the vectorization can be de done because of the automatic vectorization support of the compiler.

