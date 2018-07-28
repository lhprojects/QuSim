#define _CRT_SECURE_NO_WARNINGS
#include "SplittingMethod1DCUDA.h"

#ifdef USE_CUDA
#include "SplittingMethodGeneralCUDA.h"


template<class Scalar>
struct SplittingMethod1DCUDA : EvolverImpl1D {

	SplittingMethodGeneralCUDAImpl<Scalar> fImpl;


	virtual void initSystem1D(std::function<Complex(Real)> const & psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real)> const &vs, Real x0, Real x1, size_t nx,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, std::map<std::string, std::string> const &opts) 
	{
		EvolverImpl1D::initSystem1D(psi, force_normalization,
			dt, force_normalization_each_step,
			vs, x0, x1, nx, b, solver,
			mass, hbar, opts);

		fImpl.initSystem2D(this, fPsi.data(), force_normalization, dt, force_normalization_each_step,
			fV.data(), x0, x1, nx, 0, 1, 1,
			b, solver, mass, hbar, opts);
	}


	void step() override
	{
		fLastLastPsi = fLastPsi;
		fLastPsi = fPsi;

		//Real x = Norm2();
		update_psi();

		if (fFNES) {
			Scale(fPsi, 1.0 / sqrt(Norm2()));
		}
		fStep += fImpl.fBatch;
	}

	void update_psi() override
	{
		fImpl.update_psi();
	}


};

#endif



EvolverImpl *CreateSplittingMethod1DCUDA(std::map<std::string, std::string> const &opts)
{
#ifdef USE_CUDA
	//return new SplittingMethod2DCUDAImpl<cuComplex>();
	if (opts.find("cuda_precision") != opts.end() && opts.find("cuda_precision")->second == "single") {
		return new SplittingMethod1DCUDA<cuComplex>();
	} else {
		return new SplittingMethod1DCUDA<cuDoubleComplex>();
	}
#else
	throw std::runtime_error("cuda not supported");
#endif
}
