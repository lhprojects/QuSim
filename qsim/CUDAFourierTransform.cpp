#ifdef USE_CUDA

#include "CUDAFourierTransform.h"
#include <cuda_runtime.h>
#include <cufft.h>

#define BATCH 1

struct CUDAFourierTransform : FourierTransform {

	CUDAFourierTransform(size_t n, bool inverse);
	void Transform(std::complex<double> const *from, std::complex<double> *to) override;
	~CUDAFourierTransform();


	size_t fn;
	cufftHandle plan;
	cufftDoubleComplex *data;
	bool inv;

};

FourierTransform *CreateCUDAFourierTransform(size_t n, bool inverse)
{
	return new CUDAFourierTransform(n, inverse);
}


CUDAFourierTransform::CUDAFourierTransform(size_t n, bool inverse)
{
	inv = inverse;
	fn = n;
	int size = sizeof(cufftDoubleComplex);
	cudaError err = cudaMalloc((void**)&data, sizeof(cufftDoubleComplex)*n*BATCH);
	if (err != CUFFT_SUCCESS) {
		throw std::runtime_error("plan creation failed");
	}

	if (cudaGetLastError() != cudaSuccess) {
		throw std::runtime_error("cudaMalloc failed");
	}

	if (cufftPlan1d(&plan, (int)fn, CUFFT_Z2Z, BATCH) != CUFFT_SUCCESS) {
		throw std::runtime_error("plan creation failed");
	}

}

CUDAFourierTransform::~CUDAFourierTransform()
{
	cufftDestroy(plan);
	cudaFree(data);
}


void CUDAFourierTransform::Transform(std::complex<double> const * from, std::complex<double>* to)
{
	cudaError err = cudaMemcpy(data, from, sizeof(std::complex<double>)*fn*BATCH, cudaMemcpyKind::cudaMemcpyDefault);
	if (err != CUFFT_SUCCESS) {
		throw std::runtime_error("copy failed");
	}

	if (inv) {
		if (cufftExecZ2Z(plan, data, data, CUFFT_INVERSE) != CUFFT_SUCCESS) {
			throw std::runtime_error("ExecC2C Inverse failed");
		}
	} else {
		if (cufftExecZ2Z(plan, data, data, CUFFT_FORWARD) != CUFFT_SUCCESS) {
			throw std::runtime_error("ExecC2C Inverse failed");
		}
	}

	if (cudaDeviceSynchronize() != cudaSuccess) {
		throw std::runtime_error("Failed to synchronize");
	}
	if (cudaMemcpy(to, data, sizeof(std::complex<double>)*fn*BATCH, cudaMemcpyKind::cudaMemcpyDefault) != CUFFT_SUCCESS) {
		throw std::runtime_error("copy failed");
	}

}

#endif
