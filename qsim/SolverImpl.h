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
		std::map<std::string, std::string> const &opts);

	void MainLoop();
	void MainLoopSamllRoundError();
	void Compute() override;

	std::function<Complex(Real)> fVFunc;
	size_t fNPoints;
	size_t fNBins;
	Real fDx;
	Real fX0;
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
	Eigen::Matrix2d fTMat;

};
