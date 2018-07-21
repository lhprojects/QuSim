#pragma once
#include "EvolverImpl.h"
#include "FourierTransform.h"
#include <memory>


struct SplittingMethod : EvolverImpl1D {

	std::vector<Complex> fExpV0Dot5Dt;
	std::vector<Complex> fVPsi;
	std::vector<Complex> fTVPsi;
	std::vector<Complex> fVTVPsi;

	// period only
	std::vector<Complex> fFTPsi;
	
	std::shared_ptr<FourierTransform> fft_N;
	std::shared_ptr<FourierTransform> inv_fft_N;
	// infinite wall
	std::shared_ptr<FourierTransform> inv_fft_2N;
	std::vector<Complex> fIWPsi;
	std::vector<Complex> fIWKPsi;


	SplittingMethod();

	void initSystem1D(std::function<Complex(Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real)> const &vs, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, std::map<std::string, std::string> const &opts) override;

	void update_psi() override;
	Real CalKinEn() override;

	void initExpV();
	// vpsi = exp(-i/hbar V Dt) psi
	void ExpV(PsiVector &vpsi, PsiVector const &psi, Real t);
	void ExpT(PsiVector &tpsi, PsiVector const &psi, Real t);


};
