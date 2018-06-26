#pragma once
#include "SystemImpl.h"

struct SystemHalfVTHalfV : SystemImpl {

	std::vector<Complex> fExpV0Dot5Dt;
	std::vector<Complex> fVPsi;
	std::vector<Complex> fTVPsi;
	std::vector<Complex> fVTVPsi;

	// extend infinity only
	std::vector<Complex> fProp;
	void *inv_fft_2N;
	// period only
	std::vector<Complex> fFTPsi;
	void *fft_N;
	void *inv_fft_N;
	// infinite wall
	std::vector<Complex> fIWPsi;
	std::vector<Complex> fIWKPsi;


	SystemHalfVTHalfV()
	{
		fN = 0;
		fft_N = nullptr;
		inv_fft_N = nullptr;
		inv_fft_2N = nullptr;
	}

	void init(char const *psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		char const *vs, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar) override;

	void update_psi() override;
	Real CalKinEn() override;

	void initFreeParticleProp();
	void initExpV();
	// vpsi = exp(-i/hbar V Dt) psi
	void ExpV(PsiVector &vpsi, PsiVector const &psi, double t);
	void ExpT(PsiVector &tpsi, PsiVector const &psi);


};
