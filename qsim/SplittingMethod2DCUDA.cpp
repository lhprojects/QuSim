#define _CRT_SECURE_NO_WARNINGS
#include "EvolverImpl.h"

#ifdef USE_CUDA
#include "SplittingMethodGeneralCUDA.h"


template<class Scalar>
struct SplittingMethod2DCUDA : EvolverImpl2D {

	SplittingMethodGeneralCUDAImpl<Scalar> fImpl;


	void initSystem2D(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real, Real)> const &vs, Real x0, Real x1, size_t nx,
		Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, OptionsImpl const &opts) override
	{
		EvolverImpl2D::initSystem2D(psi, force_normalization, dt, force_normalization_each_step,
			vs, x0, x1, nx, y0, y1, ny,
			b, solver, mass, hbar, opts);
		fImpl.initSystem2D(this, fPsi.data(), force_normalization, dt, force_normalization_each_step,
			fV.data(), x0, x1, nx, y0, y1, ny,
			b, solver, mass, hbar, opts);
	}


	void step() override
	{
		fLastLastPsi = fLastPsi;
		fLastPsi = fPsi;

		//Real x = Norm2();
		update_psi();

		if (fFNES) {
			fPsi *= 1.0 / sqrt(Norm2());
		}
		fStep += fImpl.fBatch;
	}

	void update_psi() override
	{
		fImpl.update_psi();
	}

	Real CalKinEn() const override
	{
		return fImpl.CalKinEn();
	}

};

#endif



EvolverImpl2D *CreateSplittingMethod2DCUDA(OptionsImpl const &opts)
{
#ifdef USE_CUDA
	//return new SplittingMethod2DCUDAImpl<cuComplex>();
	if (opts.find("cuda_precision") != opts.end() && opts.find("cuda_precision")->second == "single") {
		return new SplittingMethod2DCUDA<cuComplex>();
	} else {
		return new SplittingMethod2DCUDA<cuDoubleComplex>();
	}
#else
	throw std::runtime_error("cuda not supported");
#endif
}
