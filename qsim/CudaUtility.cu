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

struct myf_multiplies {
	typedef cuComplex first_argument_type;
	typedef cuComplex second_argument_type;
	typedef cuComplex result_type;
	__host__ __device__ cuComplex operator()(const cuComplex &lhs, const cuComplex &rhs) const
	{
		return cuCmulf(lhs, rhs);
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

struct myf_scale {
	typedef cuComplex first_argument_type;
	typedef cuComplex second_argument_type;
	typedef cuComplex result_type;

	myf_scale(float m) : scale(m) {}
	float scale;
	__host__ __device__ cuComplex operator()(const cuComplex &lhs) const
	{
		return make_cuComplex(cuCrealf(lhs)*scale, cuCimagf(lhs)*scale);
	}
};

cudaError cudaProduct(cuDoubleComplex * __restrict__ out, const cuDoubleComplex * __restrict__ f2, const cuDoubleComplex * __restrict__ f1, size_t n)
{
	thrust::device_ptr<const cuDoubleComplex> d_f1 = thrust::device_pointer_cast(f1);
	thrust::device_ptr<const cuDoubleComplex> d_f1_e = d_f1 + n;
	thrust::device_ptr<const cuDoubleComplex> d_f2 = thrust::device_pointer_cast(f2);
	thrust::device_ptr<cuDoubleComplex> d_out = thrust::device_pointer_cast(out);

	//thrust::transform(d_f1, d_f1_e, d_f2, d_out, my_multiplies());
	return cudaSuccess;
}

cudaError cudaProduct(cuComplex * __restrict__ out, const cuComplex * __restrict__ f2, const cuComplex * __restrict__ f1, size_t n)
{
	thrust::device_ptr<const cuComplex> d_f1 = thrust::device_pointer_cast(f1);
	thrust::device_ptr<const cuComplex> d_f1_e = d_f1 + n;
	thrust::device_ptr<const cuComplex> d_f2 = thrust::device_pointer_cast(f2);
	thrust::device_ptr<cuComplex> d_out = thrust::device_pointer_cast(out);

	thrust::transform(d_f1, d_f1_e, d_f2, d_out, myf_multiplies());
	return cudaSuccess;
}


cudaError cudaScale(cuDoubleComplex * __restrict__ out, const cuDoubleComplex * __restrict__ in, double alpha, size_t n)
{
	thrust::device_ptr<const cuDoubleComplex> d_in = thrust::device_pointer_cast(in);
	thrust::device_ptr<const cuDoubleComplex> d_in_e = d_in + n;
	thrust::device_ptr<cuDoubleComplex> d_out = thrust::device_pointer_cast(out);

	thrust::transform(d_in, d_in_e, d_out, my_scale(alpha));
	return cudaSuccess;
}


cudaError cudaScale(cuComplex * __restrict__ out, const cuComplex * __restrict__ in, float alpha, size_t n)
{
	thrust::device_ptr<const cuComplex> d_in = thrust::device_pointer_cast(in);
	thrust::device_ptr<const cuComplex> d_in_e = d_in + n;
	thrust::device_ptr<cuComplex> d_out = thrust::device_pointer_cast(out);

	thrust::transform(d_in, d_in_e, d_out, myf_scale(alpha));
	return cudaSuccess;
}
