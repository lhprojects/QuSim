
#include <cuda_runtime.h>
#include <cuComplex.h>

cudaError cudaProduct(cuDoubleComplex * __restrict__ out, const cuDoubleComplex * __restrict__ f2, const cuDoubleComplex * __restrict__ f1, size_t n);
cudaError cudaScale(cuDoubleComplex * __restrict__ out, const cuDoubleComplex * __restrict__ in, double alpha, size_t n);

cudaError cudaProduct(cuComplex * __restrict__ out, const cuComplex * __restrict__ f2, const cuComplex * __restrict__ f1, size_t n);
cudaError cudaScale(cuComplex * __restrict__ out, const cuComplex * __restrict__ in, float alpha, size_t n);
