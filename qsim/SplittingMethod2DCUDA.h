#pragma once
#ifdef USE_CUDA

#include "EvolverImpl.h"
#include <cufft.h>

struct SplittingMethod2DCUDA : EvolverImpl2D {

	size_t fBatch;
	cuDoubleComplex *fExpV0Dot5Dt;
	cuDoubleComplex *fExpVC1Dt;
	cuDoubleComplex *fExpVC2Dt;

	cuDoubleComplex *fExpTDt;
	cuDoubleComplex *fExpTD1Dt;
	cuDoubleComplex *fExpTD2Dt;

	cuDoubleComplex *fTmp1;
	cuDoubleComplex *fTmp2;


	cufftHandle plan;
	Real fD1;
	Real fD2;
	Real fC1;
	Real fC2;

	SplittingMethod2DCUDA();

	void initSystem2D(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real, Real)> const &vs, Real x0, Real x1, size_t nx,
		Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, std::map<std::string, std::string> const &opts) override;

	void step() override;
	void update_psi() override;
	Real CalKinEn() const override;

	~SplittingMethod2DCUDA();

};

#endif
