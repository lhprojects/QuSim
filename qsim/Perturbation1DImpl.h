#pragma once

#include "ScatteringSolverImpl.h"
#include "FourierTransform.h"
#include "Linear.h"
#include "Perburbation.h"


struct QuPerturbation1DImpl : ScatteringSolver1DImpl {

	QuPerturbation1DImpl() : fEpsilon(), fOrder(), fSplit(), fPreconditional(), fSlow() { }

	virtual void InitPerturbation1D(
		std::function<Complex(Real)> const & v,
		Real x0,
		Real x1,
		size_t n,
		Real en,
		Real epsilon,
		Real direction,
		SolverMethod met,
		Real mass,
		Real hbar,
		std::map<std::string, std::string> const &opts);

	void Compute() override;

	Real GetMaxEnergy()
	{
		return 0.5 * pow(2 * Pi / fDx * fHbar, 2) / fMass;
	}

	Real GetMaxMomentum()
	{
		return 2 * Pi / fDx * fHbar;
	}

	Real GetMomentumGap()
	{
		return 2 * Pi / (fDx*fNx)*fHbar;
	}

	Real GetEpsilonMomentumWidth()
	{
		return fMass * abs(fEpsilon) / sqrt(2 * fMass*fE);
	}

	Real GetEnergyGap()
	{
		return sqrt(2 * fMass*fE) / fMass * 2 * Pi / (fDx*fNx)*fHbar;
	}
	Real GetEpsilonBoundaryError()
	{
		return exp(-(fNx*fDx)*sqrt(2*fMass*fE)*fEpsilon/fE);
	}


	std::shared_ptr<FourierTransform> fFFT;
	std::shared_ptr<FourierTransform> fInvFFT;

	int const fOrder;
	Real const fEpsilon;
	PsiVector fPsiK;

	int const fSplit;

	bool const fPreconditional;
	Real const fSlow;
	BornSerisePreconditioner fPreconditioner;
	std::vector<Real> fVasb;

	PsiVector ftmp1;
	PsiVector ftmp2;

};
