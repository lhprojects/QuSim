#pragma once
#include "SystemImpl.h"

// https://en.wikipedia.org/wiki/Symplectic_integrator#Splitting_methods_for_separable_Hamiltonians

struct SplittingMethod : SystemImpl1D {

	std::vector<Complex> fExpV0Dot5Dt;
	std::vector<Complex> fVPsi;
	std::vector<Complex> fTVPsi;
	std::vector<Complex> fVTVPsi;

	// period only
	std::vector<Complex> fFTPsi;
	void *fft_N;
	void *inv_fft_N;
	// infinite wall
	void *inv_fft_2N;
	std::vector<Complex> fIWPsi;
	std::vector<Complex> fIWKPsi;


	SplittingMethod()
	{
		fN = 0;
		fft_N = nullptr;
		inv_fft_N = nullptr;
		inv_fft_2N = nullptr;
	}

	void initSystem1D(char const *psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		char const *vs, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar) override;

	void update_psi() override;
	Real CalKinEn() override;

	void initExpV();
	// vpsi = exp(-i/hbar V Dt) psi
	void ExpV(PsiVector &vpsi, PsiVector const &psi, Real t);
	void ExpT(PsiVector &tpsi, PsiVector const &psi, Real t);


};
