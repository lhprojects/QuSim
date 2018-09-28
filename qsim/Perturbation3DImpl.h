#pragma once

#include "ScatteringSolverImpl.h"
#include "FourierTransform.h"
#include "FourierTransformOptions.h"
#include "Linear.h"
#include "Perburbation.h"
#include "PerturbationOptions.h"


struct QuPerturbation3DImpl : ScatteringSolver3DImpl {

	QuPerturbation3DImpl() : fEpsilon(), fOrder()
	{
	}

	virtual void InitPerturbation3D(
		std::function<Complex(Real, Real, Real)> const & v,
		Real x0,
		Real x1,
		size_t nx,
		Real y0,
		Real y1,
		size_t ny,
		Real z0,
		Real z1,
		size_t nz,
		Real en,
		Real epsilon,
		Real directionx,
		Real directiony,
		Real directionz,
		SolverMethod met,
		Real mass,
		Real hbar,
		std::map<std::string, std::string> const &opts);


	void Compute() override;

	FourierTransformOptions fFourierTransformOptions;
	std::shared_ptr<FourierTransform3D> fFFT;
	std::shared_ptr<FourierTransform3D> fInvFFT;

	PerturbationOptions fPerturbationOptions;
	int const fOrder;
	Real const fEpsilon;
	PsiVector fPsiK;


	std::vector<Real> fVasb;
	Real fNormDeltaPsi;
	PsiVector ftmp1;
	PsiVector ftmp2;

};
