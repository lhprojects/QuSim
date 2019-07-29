#include "SplittingMethod2D.h"
//#include "kissfft.hh"

void SplittingMethod2D::initSystem2D(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	std::function<Complex(Real, Real)> const &vs, Real x0, Real x1,
	size_t nx, Real y0, Real y1,
	size_t ny, BoundaryCondition b,
	SolverMethod solver, Real mass, Real hbar,
	OptionsImpl const &opts)
{
	EvolverImpl2D::initSystem2D(psi, force_normalization, dt, force_normalization_each_step,
		vs, x0, x1, nx, y0, y1, ny,
		b, solver, mass, hbar, opts);

	fVPsi.resize(fNy, fNx);
	fTVPsi.resize(fNy, fNx);
	fVTVPsi.resize(fNy, fNx);
	fPsiYIn.resize(fNx);
	fPsiYOut.resize(fNx);

	initExpV();
	initExpT();

	fFourierTransformOptions.Init(opts);

	if (b == BoundaryCondition::Period) {
		fFTPsi.resize(fNy, fNx);

		fft.reset(FourierTransform2D::Create(fNy, fNx, false, fFourierTransformOptions.fLib));
		inv_fft.reset(FourierTransform2D::Create(fNy, fNx, true, fFourierTransformOptions.fLib));

	} else {
		throw std::runtime_error("unsupported boundary condition!");
	}

	//double x = Norm2();

}

void SplittingMethod2D::update_psi()
{
	if (SolverMethod::SplittingMethodO2 == fSolverMethod) {

		ExpV(fVPsi, fPsi, 0.5);
		//double x1 = fVPsi.squaredNorm();
		ExpT(fTVPsi, fVPsi, 1.0);
		//Copy(fTVPsi, fVPsi);
		ExpV(fPsi, fTVPsi, 0.5);
		//Copy(fPsi, fVTVPsi);
		//double x2 = Norm2();
	} else if(SolverMethod::SplittingMethodO4 == fSolverMethod) {

		Real tpow1t = pow(2, 1 / 3.0);
		Real c1 = 1 / (2 * (2 - tpow1t));
		Real c2 = (1 - tpow1t) / (2 * (2 - tpow1t));

		ExpV(fVPsi, fPsi, c1);
		ExpT(fTVPsi, fVPsi, fD1);
		ExpV(fVPsi, fTVPsi, c2);
		ExpT(fTVPsi, fVPsi, fD2);
		ExpV(fVPsi, fTVPsi, c2);
		ExpT(fTVPsi, fVPsi, fD1);
		ExpV(fPsi, fTVPsi, c1);
	} else {
		throw std::runtime_error("unspported method");
	}

}

void SplittingMethod2D::initExpV()
{
	if (SolverMethod::SplittingMethodO2 == fSolverMethod) {
		fExpV0Dot5Dt.resize(fNy, fNx);

		Complex f = -1.0 / fHbar * fDt * 0.5;
		for (size_t i = 0; i < fNx*fNy; ++i) {
			fExpV0Dot5Dt.data()[i] = exp(f*fV.data()[i] * I);
		}
	}
}

void SplittingMethod2D::initExpT()
{
	if (SolverMethod::SplittingMethodO2 == fSolverMethod) {
		fExpTDt.resize(fNy, fNx);
	} else if (SolverMethod::SplittingMethodO4 == fSolverMethod) {
		fExpTD1Dt.resize(fNy, fNx);
		fExpTD2Dt.resize(fNy, fNx);
	}

	Real tpow1t = pow(2, 1 / 3.0);
	fD1 = 1 / (2 - tpow1t);
	fD2 = -tpow1t / (2 - tpow1t);

	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			size_t ii = i > fNx / 2 ? fNx - i : i;
			size_t jj = j > fNy / 2 ? fNy - j : j;
			Real kx = ii * 2 * Pi / (fDx * fNx);
			Real ky = jj * 2 * Pi / (fDy * fNy);

			Real t = fHbar * (kx*kx + ky * ky) / (2 * fMass);

			if (SolverMethod::SplittingMethodO2 == fSolverMethod) {
				fExpTDt(j, i) = exp(-I * (t * fDt* 1.0));
			} else if (SolverMethod::SplittingMethodO4 == fSolverMethod) {
				fExpTD1Dt(j, i) = exp(-I * (t * fDt* fD1));
				fExpTD2Dt(j, i) = exp(-I * (t * fDt* fD2));
			}

		}
	}

}

void SplittingMethod2D::ExpV(Eigen::MatrixXcd &vpsi, Eigen::MatrixXcd const &psi, Real t)
{
	if (t == 0.5) { // because exp() is very slow
		for (size_t i = 0; i < fNx*fNy; ++i) {
			vpsi.data()[i] = psi.data()[i] * fExpV0Dot5Dt.data()[i];
		}
	} else {
		Complex f = -1.0 / fHbar * fDt * t;
		for (size_t i = 0; i < fNx*fNy; ++i) {
			vpsi.data()[i] = psi.data()[i] * exp(f*fV.data()[i] * I);
		}
	}

}

void SplittingMethod2D::ExpT(Eigen::MatrixXcd &tpsi, Eigen::MatrixXcd const &psi, Real tt)
{
	//double y1 = psi.squaredNorm()*fDx*fDy;
	fft->Transform(psi.data(), fFTPsi.data());
	//double y2 = sqrt(fFTPsi.squaredNorm()*fDx*fDy);
	if (tt == 1.0) {
		for (size_t i = 0; i < fNx; ++i) {
			for (size_t j = 0; j < fNy; ++j) {
				fFTPsi(j, i) *= fExpTDt(j, i);
			}
		}
	} else if (tt == fD1) {
		for (size_t i = 0; i < fNx; ++i) {
			for (size_t j = 0; j < fNy; ++j) {
				fFTPsi(j, i) *= fExpTD1Dt(j, i);
			}
		}
	} else if(tt == fD2){
		for (size_t i = 0; i < fNx; ++i) {
			for (size_t j = 0; j < fNy; ++j) {
				fFTPsi(j, i) *= fExpTD2Dt(j, i);
			}
		}
	} else {
		for (size_t i = 0; i < fNx; ++i) {
			for (size_t j = 0; j < fNy; ++j) {
				size_t ii = i > fNx / 2 ? fNx - i : i;
				size_t jj = j > fNy / 2 ? fNy - j : j;
				Real kx = ii * 2 * Pi / (fDx * fNx);
				Real ky = jj * 2 * Pi / (fDy * fNy);

				Real t = fHbar * (kx*kx + ky * ky) / (2 * fMass);
				Complex f = exp(-I * (t * fDt* tt));
				fFTPsi(j, i) *= f;
			}
		}
	}


	inv_fft->Transform(fFTPsi.data(), tpsi.data());
	//double y5 = sqrt(tpsi.squaredNorm()*fDx*fDy);
	tpsi *= 1.0 / (fNx * fNy);

	//double y3 = sqrt(tpsi.squaredNorm()*fDx*fDy);

}

Real SplittingMethod2D::CalKinEn() const
{
	fft->Transform(fPsi.data(), fFTPsi.data());

	Real en = 0;

	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			size_t ii = i > fNx / 2 ? fNx - i : i;
			size_t jj = j > fNy / 2 ? fNy - j : j;
			Real kx = ii * 2 * Pi / (fDx * fNx);
			Real ky = jj * 2 * Pi / (fDy * fNy);

			Real e = fHbar*fHbar * (kx*kx + ky * ky) / (2 * fMass);
			en += abs2(fFTPsi(j, i)) * e;
		}
	}

	en /= fFTPsi.squaredNorm();
	return en;
}
