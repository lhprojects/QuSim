#pragma once

#include "SystemImpl.h"
#include "eigen/Eigen/Dense"
#include <memory>
#include "kissfft.hh"


struct SplittingMethod2D : SystemImpl2D {
	

	Eigen::MatrixXcd fExpV0Dot5Dt;
	Eigen::MatrixXcd fVPsi;
	Eigen::MatrixXcd fTVPsi;
	Eigen::MatrixXcd fVTVPsi;

	// period only
	mutable Eigen::MatrixXcd fFTPsi;
	mutable Eigen::VectorXcd fPsiYIn;
	mutable Eigen::VectorXcd fPsiYOut;
	std::shared_ptr<kissfft<Real> > fft_Nx;
	std::shared_ptr<kissfft<Real> > inv_fft_Nx;
	std::shared_ptr<kissfft<Real> > fft_Ny;
	std::shared_ptr<kissfft<Real> > inv_fft_Ny;


	SplittingMethod2D()
	{
		fNx = 0;
		fNy = 0;
	}

	void initSystem2D(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real, Real)> const &vs, Real x0, Real x1, size_t nx,
		Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, std::map<std::string, std::string> const &opts) override;

	void update_psi() override;
	Real CalKinEn() const override;

private:
	void initExpV();
	// vpsi = exp(-i/hbar V Dt) psi
	void ExpV(Eigen::MatrixXcd &vpsi, Eigen::MatrixXcd const &psi, Real t);
	void ExpT(Eigen::MatrixXcd &tpsi, Eigen::MatrixXcd const &psi, Real tt);


};
