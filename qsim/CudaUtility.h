
#include <cuda_runtime.h>
#include <cuComplex.h>

cudaError cudaProduct(cuDoubleComplex *out, cuDoubleComplex *f2, cuDoubleComplex *f1, size_t n);
cudaError cudaScale(cuDoubleComplex *out, cuDoubleComplex *in, double alpha, size_t n);
