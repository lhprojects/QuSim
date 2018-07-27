#define _CRT_SECURE_NO_WARNINGS
#include "EvolverImpl.h"

#ifdef USE_CUDA
#include "SplittingMethod2DCUDA.h"
#include "CudaUtility.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas.h>
#include <cufft.h>

template<class Scalar>
struct SplittingMethod2DCUDAImpl : EvolverImpl2D {

	Scalar *fExpVDtScaled;
	Scalar *fExpV0Dot5Dt;
	Scalar *fExpVC1Dt;
	Scalar *fExpVC2Dt;

	Scalar *fExpTDt;
	Scalar *fExpTD1Dt;
	Scalar *fExpTD2Dt;

	Scalar *fTmp1;
	Scalar *fTmp2;


	size_t fBatch;
	cufftHandle plan;
	Real fD1;
	Real fD2;
	Real fC1;
	Real fC2;

	SplittingMethod2DCUDAImpl();

	void initSystem2D(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real, Real)> const &vs, Real x0, Real x1, size_t nx,
		Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, std::map<std::string, std::string> const &opts) override;

	void step() override;
	void update_psi() override;
	Real CalKinEn() const override;

	Eigen::MatrixXcf tmpf;

	~SplittingMethod2DCUDAImpl();

};


cufftResult fftExec(cufftHandle plan,
	cufftComplex *idata,
	cufftComplex *odata,
	int direction)
{
	return cufftExecC2C(plan, idata, odata, direction);
}


cufftResult fftExec(cufftHandle plan,
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
SplittingMethod2DCUDAImpl<Scalar>::SplittingMethod2DCUDAImpl()
{
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

std::string my_itoa(int i)
{
	char b[10];
	sprintf(b, "%d", i);
	return b;
}
#define check_err(err) do { if(err!=cudaSuccess) throw std::runtime_error("cuda error: " #err " : " + my_itoa((int)err)); } while(0)


template<class Scalar>
void SplittingMethod2DCUDAImpl<Scalar>::initSystem2D(std::function<Complex(Real, Real)> const & psi,
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
		check_err(cudaMemcpy(fExpTDt, tmp1.data(), fNx*fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
	} else {
		check_err(cudaMemcpy(fExpTD1Dt, tmp1.data(), fNx*fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
		check_err(cudaMemcpy(fExpTD2Dt, tmp2.data(), fNx*fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
	}

	Complex f = -1.0 / fHbar * fDt;
	for (size_t i = 0; i < fNx*fNy; ++i) {
		if (solver  == SolverMethod::SplittingMethodO2) {
			tmp1.data()[i] = exp(f*fV.data()[i] * I*0.5);
			tmp3.data()[i] = exp(f*fV.data()[i] * I) / (1.*fNx *fNy);
		} else {
			tmp1.data()[i] = exp(f*fV.data()[i] * I*fC1);
			tmp2.data()[i] = exp(f*fV.data()[i] * I*fC2);
			tmp3.data()[i] = exp(f*fV.data()[i] * I*fC1*2.0) / (1.*fNx *fNy);
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
	check_err(cufftPlan2d(&plan, (int)fNx, (int)fNy, type));
}


template<class Scalar>
void SplittingMethod2DCUDAImpl<Scalar>::step()
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

template<class Scalar>
void SplittingMethod2DCUDAImpl<Scalar>::update_psi()
{
	if (fStep == 0) {

		if (Trait<Scalar>::Double) {
			check_err(cudaMemcpy(fTmp1, fPsi.data(), fNx * fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
		} else {
			tmpf.resize(fNy, fNx);
			tmpf = fPsi.cast<std::complex<float> >();
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
		check_err(cudaMemcpy(fPsi.data(), fTmp1, fNx * fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
	} else {
		check_err(cudaMemcpy(tmpf.data(), fTmp1, fNx * fNy * sizeof(Scalar), cudaMemcpyKind::cudaMemcpyDefault));
		fPsi = tmpf.cast<std::complex<double> >();
	}


}

template<class Scalar>
Real SplittingMethod2DCUDAImpl<Scalar>::CalKinEn() const
{
	return 0;
}

template<class Scalar>
SplittingMethod2DCUDAImpl<Scalar>::~SplittingMethod2DCUDAImpl()
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



EvolverImpl2D *CreateSplittingMethod2DCUDA(std::map<std::string, std::string> const &opts)
{
#ifdef USE_CUDA
	//return new SplittingMethod2DCUDAImpl<cuComplex>();
	if (opts.find("cuda_precision") != opts.end() && opts.find("cuda_precision")->second == "single") {
		return new SplittingMethod2DCUDAImpl<cuComplex>();
	} else {
		return new SplittingMethod2DCUDAImpl<cuDoubleComplex>();
	}
#else
	throw std::runtime_error("cuda not supported");
#endif
}
