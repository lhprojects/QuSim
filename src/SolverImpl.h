#pragma once

#include "QuSim.h"
#include "eigen/Eigen/Dense"
#include "IVPSolverImpl.h"

struct SolverImpl1D : IVPSolverImpl {

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
		Real hbar,
		OptionsImpl const &opts);

	void MainLoop();
	void MainLoopSamllRoundError();
	void Compute() override;
	void CalculateFinalJFromPsi();

	std::function<Complex(Real)> fVFunc;
	size_t const fNPoints = 0;
	size_t const fNBins = 0;
	Real const fDx = 0;
	Real const fX0 = 0;
	Real fInitJ;
	Real fFinalJ;
	PsiVector fPsi;
	PsiVector fPsiPrime;
	bool fSmallRoundError;
	std::vector<Real> fV;
	Real fV0;
	Real fV1;

	Real fT;
	Real fR;
	Eigen::MatrixXd fTMat;

};
