#pragma once

#include "ScatteringSolverImpl.h"
#include "FourierTransform.h"
#include "FourierTransformOptions.h"
#include "Linear.h"
#include "Perburbation.h"
#include "PerturbationOptions.h"


struct QuPerturbation2DImpl : ScatteringSolver2DImpl {

	QuPerturbation2DImpl() : fEpsilon(), fOrder() {
	}

	virtual void InitPerturbation2D(
		std::function<Complex(Real, Real)> const & v,
		Real x0,
		Real x1,
		size_t nx,
		Real y0,
		Real y1,
		size_t ny,
		Real en,
		Real epsilon,
		Real directionx,
		Real directiony,
		SolverMethod met,
		Real mass,
		Real hbar,
		OptionsImpl const &opts);


	void Compute() override;

	FourierTransformOptions fFourierTransformOptions;
	std::shared_ptr<FourierTransform2D> fFFT;
	std::shared_ptr<FourierTransform2D> fInvFFT;

	PerturbationOptions fPerturbationOptions;
	int const fOrder;
	Real const fEpsilon;
	PsiVector fPsiK;


	std::vector<Real> fVasb;

	PsiVector ftmp1;
	PsiVector ftmp2;

};
