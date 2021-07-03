
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>

cudaError cudaProduct(cuDoubleComplex * __restrict__ out, const cuDoubleComplex * __restrict__ f2, const cuDoubleComplex * __restrict__ f1, size_t n);
cudaError cudaScale(cuDoubleComplex * __restrict__ out, const cuDoubleComplex * __restrict__ in, double alpha, size_t n);

cudaError cudaProduct(cuComplex * __restrict__ out, const cuComplex * __restrict__ f2, const cuComplex * __restrict__ f1, size_t n);
cudaError cudaScale(cuComplex * __restrict__ out, const cuComplex * __restrict__ in, float alpha, size_t n);

inline cufftResult fftExec(cufftHandle plan,
    cufftComplex* idata,
    cufftComplex* odata,
    int direction)
{
    return cufftExecC2C(plan, idata, odata, direction);
}


inline cufftResult fftExec(cufftHandle plan,
    cufftDoubleComplex* idata,
    cufftDoubleComplex* odata,
    int direction)
{
    return cufftExecZ2Z(plan, idata, odata, direction);
}

inline cufftResult cufftExec(cufftHandle plan,
    cufftComplex* idata,
    cufftComplex* odata,
    int direction)
{
    return cufftExecC2C(plan, idata, odata, direction);
}


inline cufftResult cufftExec(cufftHandle plan,
    cufftDoubleComplex* idata,
    cufftDoubleComplex* odata,
    int direction)
{
    return cufftExecZ2Z(plan, idata, odata, direction);
}
