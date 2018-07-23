#include "CudaUtility.h"
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/device_ptr.h>

struct my_multiplies {
	typedef cuDoubleComplex first_argument_type;
	typedef cuDoubleComplex second_argument_type;
	typedef cuDoubleComplex result_type;
	__host__ __device__ cuDoubleComplex operator()(const cuDoubleComplex &lhs, const cuDoubleComplex &rhs) const
	{
		return cuCmul(lhs, rhs);
	}
};

struct my_scale {
	typedef cuDoubleComplex first_argument_type;
	typedef cuDoubleComplex second_argument_type;
	typedef cuDoubleComplex result_type;

	my_scale(double m) : scale(m) { }
	double scale;
	__host__ __device__ cuDoubleComplex operator()(const cuDoubleComplex &lhs) const
	{
		return make_cuDoubleComplex(cuCreal(lhs)*scale, cuCimag(lhs)*scale);
	}
};

#if 0
__global__
void mul(size_t n, cuDoubleComplex *out, cuDoubleComplex *x, cuDoubleComplex *y)
{
	size_t i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i < n) out[i] = cuCmul(x[i], y[i]);
}

cudaError cudaProduct(cuDoubleComplex *out, cuDoubleComplex * f2, cuDoubleComplex * f1, size_t n)
{
	saxpy <<< (n + 255) / 256, 256 >>>(n, out, f2, f1);

	return cudaSuccess;
}

#else
cudaError cudaProduct(cuDoubleComplex *out, cuDoubleComplex * f2, cuDoubleComplex * f1, size_t n)
{
	thrust::device_ptr<cuDoubleComplex> d_f1 = thrust::device_pointer_cast(f1);
	thrust::device_ptr<cuDoubleComplex> d_f1_e = d_f1 + n;
	thrust::device_ptr<cuDoubleComplex> d_f2 = thrust::device_pointer_cast(f2);
	thrust::device_ptr<cuDoubleComplex> d_out = thrust::device_pointer_cast(out);

	thrust::transform(d_f1, d_f1_e, d_f2, d_out, my_multiplies());
	return cudaSuccess;
}
#endif


cudaError cudaScale(cuDoubleComplex *out, cuDoubleComplex *in, double alpha, size_t n)
{
	thrust::device_ptr<cuDoubleComplex> d_in = thrust::device_pointer_cast(in);
	thrust::device_ptr<cuDoubleComplex> d_in_e = d_in + n;
	thrust::device_ptr<cuDoubleComplex> d_out = thrust::device_pointer_cast(out);

	thrust::transform(d_in, d_in_e, d_out, my_scale(alpha));
	return cudaSuccess;
}
