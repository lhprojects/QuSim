#pragma once

#include "PerturbationImpl.h"
#include "FourierTransform.h"

struct QuPerturbation1DImpl : QuPerturbationImpl {

	QuPerturbation1DImpl() : fDx(), fX0(), fEpsilon(), fNx(), fK0() { }
	virtual void InitQuPerturbation1D(
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

	virtual void Compute();
	Real GetX(ptrdiff_t i) const { return fX0 + fDx * i; }
	void initPotential();

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


	Real const fEpsilon;
	size_t const fNx;
	std::function<Complex(Real)> const fVFunc;
	Real const fX0;
	Real const fDx;
	Real const fK0;

	PsiVector fPsi0X;
	PsiVector fPsi1X;
	PsiVector fPsi2X;
	PsiVector fPsi0K;
	PsiVector fPsi1K;
	PsiVector fPsi2K;

	PsiVector ftmp1;
	PsiVector ftmp2;

	PsiVector fPsi;
	std::vector<Real> fV;
	Real fT;
	Real fR;

};
