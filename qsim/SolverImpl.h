#pragma once

#include "QuSim.h"
#include "eigen/Eigen/Dense"

struct SolverImpl {

	SolverImpl();
	virtual ~SolverImpl() { }

	void initSystem(
		Real en,
		Real mass,
		Real hbar);

	virtual void Calculate() = 0;

	Real const fE;
	Real const fMass;
	Real const fHbar;

};

struct SolverImpl1D : SolverImpl {

	SolverImpl1D();

	void initSystem1D(
		std::function<Complex(Real)> const & v,
		Real x0,
		Real x1,
		size_t n,
		Real en,
		Complex initPsi,
		Complex initPsiPrime,
		SolverMethod met,
		Real mass,
		Real hbar);

	void Calculate() override;

	std::function<Complex(Real)> fVFunc;
	size_t fNPoints;
	size_t fNBins;
	Real fDx;
	Real fX0;
	Real fE;
	Real fInitJ;
	Real fFinalJ;
	PsiVector fPsi;
	PsiVector fPsiPrime;
	SolverMethod fMethod;
	std::vector<Real> fV;
	Real fMass;
	Real fHbar;
	Real fV0;
	Real fV1;

	Real fT;
	Real fR;
	Eigen::Matrix2cd fTMat;
	Real T11;
	Real T12;
	Real T21;
	Real T22;

};
