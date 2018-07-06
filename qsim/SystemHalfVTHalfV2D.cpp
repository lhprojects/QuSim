#include "SystemHalfVTHalfV2D.h"
//#include "kissfft.hh"

void SystemHalfVTHalfV2D::initSystem2D(char const * psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	char const * vs, Real x0, Real x1,
	size_t nx, Real y0, Real y1,
	size_t ny, BoundaryCondition b,
	SolverMethod solver, Real mass, Real hbar)
{
	SystemImpl2D::initSystem2D(psi, force_normalization, dt, force_normalization_each_step,
		vs, x0, x1, nx, y0, y1, ny,
		b, solver, mass, hbar);

	fVPsi.resize(fNy, fNx);
	fTVPsi.resize(fNy, fNx);
	fVTVPsi.resize(fNy, fNx);
	fPsiYIn.resize(fNx);
	fPsiYOut.resize(fNx);

	initExpV();
	if (b == BoundaryCondition::Period) {
		fFTPsi.resize(fNy, fNx);

		fft_Nx.reset(new kissfft<Real>(fNx, false));
		inv_fft_Nx.reset(new kissfft<Real>(fNx, true));
		if (fNx == fNy) {
			fft_Ny = fft_Nx;
			inv_fft_Ny = inv_fft_Nx;
		} else {
			fft_Ny.reset(new kissfft<Real>(fNy, false));
			inv_fft_Ny.reset(new kissfft<Real>(fNy, true));
		}

	} else {
		throw std::runtime_error("unsupported boundary condition!");
	}

	//double x = Norm2();

}

void SystemHalfVTHalfV2D::update_psi()
{
	ExpV(fVPsi, fPsi, 0.5);
	//double x1 = fVPsi.squaredNorm();
	ExpT(fTVPsi, fVPsi);
	//Copy(fTVPsi, fVPsi);
	ExpV(fPsi, fTVPsi, 0.5);
	//Copy(fPsi, fVTVPsi);
	//double x2 = Norm2();

}

void SystemHalfVTHalfV2D::initExpV()
{
	fExpV0Dot5Dt.resize(fNy, fNx);

	Complex f = -1.0 / fHbar * fDt * 0.5;
	for (size_t i = 0; i < fNx*fNy; ++i) {
		fExpV0Dot5Dt.data()[i] = exp(f*fV.data()[i] * I);
	}

}


void SystemHalfVTHalfV2D::ExpV(Eigen::MatrixXcd &vpsi, Eigen::MatrixXcd const &psi, double t)
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

void SystemHalfVTHalfV2D::ExpT(Eigen::MatrixXcd &tpsi, Eigen::MatrixXcd const &psi)
{
	//double x = Norm2();

	for (size_t i = 0; i < fNx; ++i) {
		fft_Ny->transform(psi.data() + i * fNy, fFTPsi.data() + i * fNy);
	}

	//double y = fFTPsi.squaredNorm()*fDx*fDy;

	for (size_t i = 0; i < fNy; ++i) {
		for (size_t j = 0; j < fNx; ++j) {
			fPsiYIn.data()[j] = fFTPsi.data()[j*fNy + i];
		}
		fft_Nx->transform(fPsiYIn.data(), fPsiYOut.data());
		for (size_t j = 0; j < fNx; ++j) {
			fFTPsi.data()[j*fNy + i] = fPsiYOut[j];
		}
	}
	
	//double y2 = fFTPsi.squaredNorm()*fDx*fDy;
	fFTPsi *= 1.0 / sqrt(fNx * fNy);


	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			size_t ii = i > fNx / 2 ? fNx - i : i;
			size_t jj = j > fNy / 2 ? fNy - j : j;
			Real kx = ii * 2 * Pi / (fDx * fNx);
			Real ky = jj * 2 * Pi / (fDy * fNy);

			Real t = fHbar * (kx*kx + ky * ky) / (2 * fMass);
			Complex f = exp(-I * (t * fDt));
			fFTPsi(j, i) *= f;
		}
	}

	for (size_t i = 0; i < fNy; ++i) {
		for (size_t j = 0; j < fNx; ++j) {
			fPsiYIn[j] = fFTPsi.data()[j*fNy + i];
		}
		inv_fft_Nx->transform(fPsiYIn.data(), fPsiYOut.data());
		for (size_t j = 0; j < fNx; ++j) {
			fFTPsi.data()[j*fNy + i] = fPsiYOut[j];
		}
	}
	for (size_t i = 0; i < fNx; ++i) {
		inv_fft_Ny->transform(fFTPsi.data() + i * fNy, tpsi.data() + i * fNy);
	}

	tpsi *= 1.0 / sqrt(fNx * fNy);

	//double x = Norm2();

}

Real SystemHalfVTHalfV2D::CalKinEn() const
{
	for (size_t i = 0; i < fNx; ++i) {
		fft_Ny->transform(fPsi.data() + i * fNy, fFTPsi.data() + i * fNy);
	}

	for (size_t i = 0; i < fNy; ++i) {
		for (size_t j = 0; j < fNx; ++j) {
			fPsiYIn[j] = fFTPsi.data()[j*fNy + i];
		}
		fft_Nx->transform(fPsiYIn.data(), fPsiYOut.data());
		for (size_t j = 0; j < fNx; ++j) {
			fFTPsi.data()[j*fNy + i] = fPsiYOut[j];
		}
	}

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