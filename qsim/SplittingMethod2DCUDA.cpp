#ifdef USE_CUDA
#define _CRT_SECURE_NO_WARNINGS
#include "SplittingMethod2DCUDA.h"
#include "CudaUtility.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas.h>

SplittingMethod2DCUDA::SplittingMethod2DCUDA()
{

	fExpV0Dot5Dt = nullptr;
	fExpVC1Dt = nullptr;
	fExpVC2Dt = nullptr;

	fExpTDt = nullptr;
	fExpTD1Dt = nullptr;
	fExpTD2Dt = nullptr;

	fTmp1 = nullptr;
	fTmp2 = nullptr;

}

#define check_err(err) do { if(err!=cudaSuccess) throw std::runtime_error("cuda error: " #err); } while(0)

void SplittingMethod2DCUDA::initSystem2D(std::function<Complex(Real, Real)> const & psi,
	bool force_normalization, Complex dt, bool force_normalization_each_step,
	std::function<Complex(Real, Real)> const & vs,
	Real x0, Real x1, size_t nx, Real y0, Real y1, size_t ny,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar, std::map<std::string, std::string> const & opts)
{
	EvolverImpl2D::initSystem2D(psi, force_normalization, dt, force_normalization_each_step,
		vs, x0, x1, nx, y0, y1, ny,
		b, solver, mass, hbar, opts);

	fBatch = 1;
	{
		auto it = opts.find("batch");
		if (it != opts.end()) {
			auto x = it->second;
			int d;
			if (sscanf(x.c_str(), "%d", &d) < 0) {
				throw std::runtime_error("batch parse error");
			}
			fBatch = d;
		}
	}

	check_err(cudaMalloc((void**)&fTmp1, nx*ny * sizeof(cuDoubleComplex)));
	check_err(cudaMalloc((void**)&fTmp2, nx*ny * sizeof(cuDoubleComplex)));

	if (solver == SolverMethod::SplittingMethodO2) {
		check_err(cudaMalloc((void**)&fExpV0Dot5Dt, nx*ny * sizeof(cuDoubleComplex)));
		check_err(cudaMalloc((void**)&fExpTDt, nx*ny * sizeof(cuDoubleComplex)));
	} else {
		check_err(cudaMalloc((void**)&fExpVC1Dt, nx*ny * sizeof(cuDoubleComplex)));
		check_err(cudaMalloc((void**)&fExpVC2Dt, nx*ny * sizeof(cuDoubleComplex)));
		check_err(cudaMalloc((void**)&fExpTD1Dt, nx*ny * sizeof(cuDoubleComplex)));
		check_err(cudaMalloc((void**)&fExpTD2Dt, nx*ny * sizeof(cuDoubleComplex)));

		Real tpow1t = pow(2, 1 / 3.0);
		fD1 = 1 / (2 - tpow1t);
		fD2 = -tpow1t / (2 - tpow1t);
		fC1 = 1 / (2 * (2 - tpow1t));
		fC2 = (1 - tpow1t) / (2 * (2 - tpow1t));

	}

	Eigen::MatrixXcd tmp1(fNy, fNx);
	Eigen::MatrixXcd tmp2(fNy, fNx);
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			size_t ii = i > fNx / 2 ? fNx - i : i;
			size_t jj = j > fNy / 2 ? fNy - j : j;
			Real kx = ii * 2 * Pi / (fDx * fNx);
			Real ky = jj * 2 * Pi / (fDy * fNy);

			Real t = fHbar * (kx*kx + ky * ky) / (2 * fMass);

			if (solver == SolverMethod::SplittingMethodO2) {
				tmp1(j, i) = exp(-I * (t * fDt));
			} else {
				tmp1(j, i) = exp(-I * (t * fDt* fD1));
				tmp2(j, i) = exp(-I * (t * fDt* fD2));
			}
		}
	}

	if (solver == SolverMethod::SplittingMethodO2) {
		check_err(cudaMemcpy(fExpTDt, tmp1.data(), fNx*fNy * sizeof(cuDoubleComplex), cudaMemcpyKind::cudaMemcpyDefault));
	} else {
		check_err(cudaMemcpy(fExpTD1Dt, tmp1.data(), fNx*fNy * sizeof(cuDoubleComplex), cudaMemcpyKind::cudaMemcpyDefault));
		check_err(cudaMemcpy(fExpTD2Dt, tmp2.data(), fNx*fNy * sizeof(cuDoubleComplex), cudaMemcpyKind::cudaMemcpyDefault));
	}

	Complex f = -1.0 / fHbar * fDt;
	for (size_t i = 0; i < fNx*fNy; ++i) {
		if (solver  == SolverMethod::SplittingMethodO2) {
			tmp1.data()[i] = exp(f*fV.data()[i] * I*0.5);
		} else {
			tmp1.data()[i] = exp(f*fV.data()[i] * I*fC1);
			tmp1.data()[i] = exp(f*fV.data()[i] * I*fC2);
		}
	}

	if (solver == SolverMethod::SplittingMethodO2) {
		check_err(cudaMemcpy(fExpV0Dot5Dt, tmp1.data(), fNx*fNy * sizeof(cuDoubleComplex), cudaMemcpyKind::cudaMemcpyDefault));
	} else {
		check_err(cudaMemcpy(fExpVC1Dt, tmp1.data(), fNx*fNy * sizeof(cuDoubleComplex), cudaMemcpyKind::cudaMemcpyDefault));
		check_err(cudaMemcpy(fExpVC2Dt, tmp2.data(), fNx*fNy * sizeof(cuDoubleComplex), cudaMemcpyKind::cudaMemcpyDefault));
	}

	check_err(cufftPlan2d(&plan, (int)fNx, (int)fNy, CUFFT_Z2Z));
}

void SplittingMethod2DCUDA::step()
{
	//Real y = Norm2();
	fLastLastPsi = fLastPsi;
	fLastPsi = fPsi;

	//Real x = Norm2();
	update_psi();

	if (fFNES) {
		fPsi *= 1.0 / sqrt(Norm2());
	}
	fStep += fBatch;

}

void SplittingMethod2DCUDA::update_psi()
{
	if (fStep == 0) {
		check_err(cudaMemcpy(fTmp1, fPsi.data(), fNx * fNy * sizeof(cuDoubleComplex), cudaMemcpyKind::cudaMemcpyDefault));
	}
	for (size_t i = 0; i < fBatch; ++i) {
		if (SolverMethod::SplittingMethodO2 == fSolverMethod) {

			check_err(cudaProduct(fTmp2, fExpV0Dot5Dt, fTmp1, fNx*fNy));
			check_err(cufftExecZ2Z(plan, fTmp2, fTmp1, CUFFT_FORWARD));
			check_err(cudaProduct(fTmp2, fExpTDt, fTmp1, fNx*fNy));
			check_err(cufftExecZ2Z(plan, fTmp2, fTmp1, CUFFT_INVERSE));
			check_err(cudaProduct(fTmp2, fExpV0Dot5Dt, fTmp1, fNx*fNy));
			check_err(cudaScale(fTmp1, fTmp2, 1. / (fNx*fNy), fNx * fNy));



		} else {

			check_err(cudaProduct(fTmp2, fExpVC1Dt, fTmp1, fNx*fNy));

			check_err(cufftExecZ2Z(plan, fTmp2, fTmp1, CUFFT_FORWARD));
			check_err(cudaProduct(fTmp2, fExpTD1Dt, fTmp1, fNx*fNy));
			check_err(cufftExecZ2Z(plan, fTmp2, fTmp1, CUFFT_INVERSE));

			check_err(cudaProduct(fTmp2, fExpVC2Dt, fTmp1, fNx*fNy));

			check_err(cufftExecZ2Z(plan, fTmp2, fTmp1, CUFFT_FORWARD));
			check_err(cudaProduct(fTmp2, fExpTD2Dt, fTmp1, fNx*fNy));
			check_err(cufftExecZ2Z(plan, fTmp2, fTmp1, CUFFT_INVERSE));

			check_err(cudaProduct(fTmp2, fExpVC2Dt, fTmp1, fNx*fNy));

			check_err(cufftExecZ2Z(plan, fTmp2, fTmp1, CUFFT_FORWARD));
			check_err(cudaProduct(fTmp2, fExpTD1Dt, fTmp1, fNx*fNy));
			check_err(cufftExecZ2Z(plan, fTmp2, fTmp1, CUFFT_INVERSE));

			check_err(cudaProduct(fTmp2, fExpVC1Dt, fTmp1, fNx*fNy));
			check_err(cudaScale(fTmp1, fTmp2, 1. / (fNx*fNy) / (fNx*fNy) / (fNx*fNy), fNx * fNy));

		}
	}

	check_err(cudaMemcpy(fPsi.data(), fTmp1, fNx * fNy * sizeof(cuDoubleComplex), cudaMemcpyKind::cudaMemcpyDefault));

}

Real SplittingMethod2DCUDA::CalKinEn() const
{
	return 0;
}

SplittingMethod2DCUDA::~SplittingMethod2DCUDA()
{
	if (fSolverMethod == SolverMethod::SplittingMethodO2) {
		cudaFree(fExpV0Dot5Dt);
		cudaFree(fExpTDt);
	} else {
		cudaFree(fExpVC1Dt);
		cudaFree(fExpVC2Dt);
		cudaFree(fExpTD1Dt);
		cudaFree(fExpTD2Dt);
	}
	cudaFree(fTmp1);
	cudaFree(fTmp2);

}

#endif
