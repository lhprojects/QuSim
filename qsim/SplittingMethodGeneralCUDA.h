#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include "EvolverImpl.h"

#ifdef USE_CUDA
#include "CudaUtility.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas.h>
#include <cufft.h>

template<class Scalar>
struct SplittingMethodGeneralCUDAImpl {

	EvolverImpl *fEvolver;

	Complex *fPsi;
	size_t fNx;
	size_t fNy;
	SolverMethod fSolverMethod;


	Scalar *fExpVDtScaled;
	Scalar *fExpV0Dot5Dt;
	Scalar *fExpVC1Dt;
	Scalar *fExpVC2Dt;

	Scalar *fExpTDt;
	Scalar *fExpTD1Dt;
	Scalar *fExpTD2Dt;

	Scalar *fTmp1;
	Scalar *fTmp2;
	Eigen::MatrixXcf tmpf;


	size_t fBatch;
	cufftHandle plan;
	Real fD1;
	Real fD2;
	Real fC1;
	Real fC2;

	SplittingMethodGeneralCUDAImpl();

	void initSystem2D(EvolverImpl * evolver2d, Complex *psi,
		bool force_normalization, Complex dt, bool force_normalization_each_step,
		Real const *fV,
		Real x0, Real x1, size_t nx, Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, std::map<std::string, std::string> const & opts);

	void update_psi();
	Real CalKinEn() const;

	~SplittingMethodGeneralCUDAImpl();

};

template<class Scalar>
struct SplittingMethod2DCUDAImpl : Evolver2D {

};

inline cufftResult fftExec(cufftHandle plan,
	cufftComplex *idata,
	cufftComplex *odata,
	int direction)
{
	return cufftExecC2C(plan, idata, odata, direction);
}


inline cufftResult fftExec(cufftHandle plan,
	cufftDoubleComplex *idata,
	cufftDoubleComplex *odata,
	int direction)
{
	return cufftExecZ2Z(plan, idata, odata, direction);
}

template<class Scalar>
struct Trait {

};

template<>
struct Trait<cuDoubleComplex> {
	typedef Eigen::MatrixXcd Matrix;
	typedef double Real;
	static const bool Double = true;
};

template<>
struct Trait<cuComplex> {
	typedef Eigen::MatrixXcf Matrix;
	typedef float Real;
	static const bool Double = false;

};

template<class Scalar>
SplittingMethodGeneralCUDAImpl<Scalar>::SplittingMethodGeneralCUDAImpl()
{
	fEvolver = nullptr;

	fExpVDtScaled = nullptr;
	fExpV0Dot5Dt = nullptr;
	fExpVC1Dt = nullptr;
	fExpVC2Dt = nullptr;

	fExpTDt = nullptr;
	fExpTD1Dt = nullptr;
	fExpTD2Dt = nullptr;

	fTmp1 = nullptr;
	fTmp2 = nullptr;

}

inline std::string my_itoa(int i)
{
	char b[10];
	sprintf(b, "%d", i);
	return b;
}
#define check_err(err) do { if(err!=cudaSuccess) throw std::runtime_error("cuda error: " #err " : " + my_itoa((int)err)); } while(0)


template<class Scalar>
void SplittingMethodGeneralCUDAImpl<Scalar>::initSystem2D(EvolverImpl *evolver2d, Complex *psi,
	bool force_normalization, Complex dt, bool force_normalization_each_step,
	Real const *fV,
	Real x0, Real x1, size_t nx, Real y0, Real y1, size_t ny,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar, std::map<std::string, std::string> const & opts)
{
	fEvolver = evolver2d;
	fNx = nx;
	fNy = ny;
	fSolverMethod = solver;
	fPsi = psi;

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

	check_err(cudaMalloc((void**)&fTmp1, nx*ny * sizeof(Scalar)));
	check_err(cudaMalloc((void**)&fTmp2, nx*ny * sizeof(Scalar)));
	check_err(cudaMalloc((void**)&fExpVDtScaled, nx*ny * sizeof(Scalar)));

	if (solver == SolverMethod::SplittingMethodO2) {
		check_err(cudaMalloc((void**)&fExpV0Dot5Dt, nx*ny * sizeof(Scalar)));
		check_err(cudaMalloc((void**)&fExpTDt, nx*ny * sizeof(Scalar)));
	} else {
		check_err(cudaMalloc((void**)&fExpVC1Dt, nx*ny * sizeof(Scalar)));
		check_err(cudaMalloc((void**)&fExpVC2Dt, nx*ny * sizeof(Scalar)));
		check_err(cudaMalloc((void**)&fExpTD1Dt, nx*ny * sizeof(Scalar)));
		check_err(cudaMalloc((void**)&fExpTD2Dt, nx*ny * sizeof(Scalar)));

		Real tpow1t = pow(2, 1 / 3.0);
		fD1 = 1 / (2 - tpow1t);
		fD2 = -tpow1t / (2 - tpow1t);
		fC1 = 1 / (2 * (2 - tpow1t));
		fC2 = (1 - tpow1t) / (2 * (2 - tpow1t));

	}

	typename Trait<Scalar>::Matrix tmp1(fNy, fNx);
	typename Trait<Scalar>::Matrix tmp2(fNy, fNx);
	typename Trait<Scalar>::Matrix tmp3(fNy, fNx);
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			size_t ii = i > fNx / 2 ? fNx - i : i;
			size_t jj = j > fNy / 2 ? fNy - j : j;
			Real kx = ii * 2 * Pi / (x1 - x0);
			Real ky = jj * 2 * Pi / (y1 - y0);

			Real t = hbar * (kx*kx + ky * ky) / (2 * mass);

			if (solver == SolverMethod::SplittingMethodO2) {
				tmp1(j, i) = exp(-I * (t * dt));
			} else {
				tmp1(j, i) = exp(-I * (t * dt* fD1));
				tmp2(j, i) = exp(-I * (t * dt* fD2));
			}
		}
	}

	if (solver == SolverMethod::SplittingMethodO2) {
		check_err(cudaMemcpy(fExpTDt, tmp1.data(), fNx*fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
	} else {
		check_err(cudaMemcpy(fExpTD1Dt, tmp1.data(), fNx*fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
		check_err(cudaMemcpy(fExpTD2Dt, tmp2.data(), fNx*fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
	}

	Complex f = -1.0 / hbar * dt;
	for (size_t i = 0; i < fNx*fNy; ++i) {
		if (solver == SolverMethod::SplittingMethodO2) {
			tmp1.data()[i] = exp(f*fV[i] * I*0.5);
			tmp3.data()[i] = exp(f*fV[i] * I) / (1.*fNx *fNy);
		} else {
			tmp1.data()[i] = exp(f*fV[i] * I*fC1);
			tmp2.data()[i] = exp(f*fV[i] * I*fC2);
			tmp3.data()[i] = exp(f*fV[i] * I*fC1*2.0) / (1.*fNx *fNy);
		}
	}

	check_err(cudaMemcpy(fExpVDtScaled, tmp3.data(), fNx*fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
	if (solver == SolverMethod::SplittingMethodO2) {
		check_err(cudaMemcpy(fExpV0Dot5Dt, tmp1.data(), fNx*fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
	} else {
		check_err(cudaMemcpy(fExpVC1Dt, tmp1.data(), fNx*fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
		check_err(cudaMemcpy(fExpVC2Dt, tmp2.data(), fNx*fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
	}

	cufftType_t type = Trait<Scalar>::Double ? CUFFT_Z2Z : CUFFT_C2C;
	if (fNy == 1) {
		check_err(cufftPlan1d(&plan, (int)fNx, type, 1));
	} else {
		check_err(cufftPlan2d(&plan, (int)fNx, (int)fNy, type));
	}
}

template<class Scalar>
void SplittingMethodGeneralCUDAImpl<Scalar>::update_psi()
{
	if (fEvolver->fStep == 0) {

		if (Trait<Scalar>::Double) {
			check_err(cudaMemcpy(fTmp1, fPsi, fNx * fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
		} else {
			tmpf.resize(fNy, fNx);
			for (size_t i = 0; i < fNx*fNy; ++i) { tmpf.data()[i] = (std::complex<float>)fPsi[i]; }
			check_err(cudaMemcpy(fTmp1, tmpf.data(), fNx * fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
		}
	}

	for (size_t i = 0; i < fBatch; ++i) {
		bool enable_head_tail_fold = true;
		if (SolverMethod::SplittingMethodO2 == fSolverMethod) {

			if (i == 0 || !enable_head_tail_fold) {
				check_err(cudaProduct(fTmp2, fExpV0Dot5Dt, fTmp1, fNx*fNy));
			} else {
				std::swap(fTmp1, fTmp2);
			}

			check_err(fftExec(plan, fTmp2, fTmp1, CUFFT_FORWARD));
			check_err(cudaProduct(fTmp2, fExpTDt, fTmp1, fNx*fNy));
			check_err(fftExec(plan, fTmp2, fTmp1, CUFFT_INVERSE));

			if (i == fBatch - 1 || !enable_head_tail_fold) {
				check_err(cudaProduct(fTmp2, fExpV0Dot5Dt, fTmp1, fNx*fNy));
				check_err(cudaScale(fTmp1, fTmp2, (typename Trait<Scalar>::Real)(1. / (fNx*fNy)), fNx * fNy));
			} else {
				check_err(cudaProduct(fTmp2, fExpVDtScaled, fTmp1, fNx*fNy));
				std::swap(fTmp1, fTmp2);
			}


		} else {

			check_err(cudaProduct(fTmp2, fExpVC1Dt, fTmp1, fNx*fNy));

			check_err(fftExec(plan, fTmp2, fTmp1, CUFFT_FORWARD));
			check_err(cudaProduct(fTmp2, fExpTD1Dt, fTmp1, fNx*fNy));
			check_err(fftExec(plan, fTmp2, fTmp1, CUFFT_INVERSE));

			check_err(cudaProduct(fTmp2, fExpVC2Dt, fTmp1, fNx*fNy));

			check_err(fftExec(plan, fTmp2, fTmp1, CUFFT_FORWARD));
			check_err(cudaProduct(fTmp2, fExpTD2Dt, fTmp1, fNx*fNy));
			check_err(fftExec(plan, fTmp2, fTmp1, CUFFT_INVERSE));

			check_err(cudaProduct(fTmp2, fExpVC2Dt, fTmp1, fNx*fNy));

			check_err(fftExec(plan, fTmp2, fTmp1, CUFFT_FORWARD));
			check_err(cudaProduct(fTmp2, fExpTD1Dt, fTmp1, fNx*fNy));
			check_err(fftExec(plan, fTmp2, fTmp1, CUFFT_INVERSE));

			check_err(cudaProduct(fTmp2, fExpVC1Dt, fTmp1, fNx*fNy));
			check_err(cudaScale(fTmp1, fTmp2, (typename Trait<Scalar>::Real)(1. / (fNx*fNy) / (fNx*fNy) / (fNx*fNy)), fNx * fNy));

		}
	}

	if (Trait<Scalar>::Double) {
		check_err(cudaMemcpy(fPsi, fTmp1, fNx * fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
	} else {
		check_err(cudaMemcpy(tmpf.data(), fTmp1, fNx * fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
		for (size_t i = 0; i < fNx*fNy; ++i) { fPsi[i] = (std::complex<double>)tmpf.data()[i]; }
	}


}

template<class Scalar>
Real SplittingMethodGeneralCUDAImpl<Scalar>::CalKinEn() const
{
	return 0;
}

template<class Scalar>
SplittingMethodGeneralCUDAImpl<Scalar>::~SplittingMethodGeneralCUDAImpl()
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
	cudaFree(fExpVDtScaled);
	cudaFree(fTmp1);
	cudaFree(fTmp2);

}

#endif

