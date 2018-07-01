
#include "SystemHalfVTHalfV.h"
#include "kissfft.hh"

void SystemHalfVTHalfV::init(char const *psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	char const *vs, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar)
{
	SystemImpl::init(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar);

	fVPsi.resize(n);
	fTVPsi.resize(n);
	fVTVPsi.resize(n);

	initExpV();
	if (fBoundaryCondition == BoundaryCondition::Period) {
		fFTPsi.resize(n);

		if (fft_N) delete (kissfft<Real>*)fft_N;
		fft_N = new kissfft<Real>(fN, false);
		if (inv_fft_N) delete (kissfft<Real>*)(inv_fft_N);
		inv_fft_N = new kissfft<Real>(fN, true);

	} else if (fBoundaryCondition == BoundaryCondition::InfiniteWall) {
		if (inv_fft_2N) delete (kissfft<Real>*)fft_N;
		inv_fft_2N = new kissfft<Real>(2 * fN, true);
		fIWPsi.resize(2 * n);
		fIWKPsi.resize(2 * n);
	}
}

void SystemHalfVTHalfV::initExpV()
{
	fExpV0Dot5Dt.resize(fN);

	Complex f = -1.0 / fHbar * fDt * 0.5;
	for (size_t i = 0; i < fN; ++i) {
		fExpV0Dot5Dt[i] = exp(f*fV[i] * I);
	}

}



void SystemHalfVTHalfV::ExpT(PsiVector &tpsi, PsiVector const &psi)
{
	Zero(tpsi);
	if (fBoundaryCondition == BoundaryCondition::Period) {
		//Zero(fFTPsi);
		//fKinEnergy = 0;

		kissfft<Real> &fft = *(kissfft<Real>*)this->fft_N;
		kissfft<Real> &inv_fft = *(kissfft<Real>*)this->inv_fft_N;
		//kiss_fft_cfg kiss_cfg =  kiss_fft_alloc(fN, false, NULL, NULL);
		//kiss_fft_cfg kiss_cfg_inv = kiss_fft_alloc(fN, true, NULL, NULL);

		fft.transform(psi.data(), fFTPsi.data());
		//kiss_fft(kiss_cfg, (kiss_fft_cpx*)psi.data(), (kiss_fft_cpx*)fFTPsi.data());

		Scale(fFTPsi, 1.0 / sqrt(1.0*fN));


		// Dk = 2 Pi / (a * N)

		// T = hbar^2 / (2 m)  (2 - 2cos(k a)) / a^2
		// T = hbar^2 / m  2 sin^2(0.5 k a) / a^2
		// T ~ hbar^2 k^2 / (2m)
		// k = i Dk
		for (size_t i = 0; i < fN; ++i) {
			size_t j = i > fN / 2 ? fN - i : i;
			Real k = j * 2 * Pi / (fDx * fN);
#if 0
			Real s = sin(0.5 * k * fDx);
			Real t = fHbar * (2 * s*s / (fDx*fDx)) / fMass;
#else
			Real t = fHbar * k*k / (2 * fMass);
#endif
			Complex f = exp(-I * (t * fDt));
			fFTPsi[i] *= f;

			//fKinEnergy += t*fHbar*abs2(fFTPsi[i]);
		}
		//fKinEnergy /= Norm2(fFTPsi);

		inv_fft.transform(fFTPsi.data(), tpsi.data());
		//kiss_fft(kiss_cfg_inv, (kiss_fft_cpx*)fFTPsi.data(), (kiss_fft_cpx*)tpsi.data());

		Scale(tpsi, 1.0 / sqrt(1.0*fN));

		//free(kiss_cfg);
		//free(kiss_cfg_inv);


	} else if (fBoundaryCondition == BoundaryCondition::InfiniteWall) {
		//fKinEnergy = 0;
		//Real kinEnergyNorm2 = 0;

		kissfft<Real> &inv_fft = *(kissfft<Real>*)this->inv_fft_2N;

		//kiss_fft_cfg kiss_cfg_inv = kiss_fft_alloc(2*fN, true, NULL, NULL);

		std::copy(psi.begin(), psi.end(), fIWPsi.begin());
		fIWPsi[0] = 0;
		fIWPsi[fN] = 0;
		for (size_t i = fN + 1; i < 2 * fN; ++i) {
			fIWPsi[i] = -fIWPsi[2 * fN - i];
		}

		inv_fft.transform(fIWPsi.data(), fIWKPsi.data());
		//kiss_fft(kiss_cfg_inv, (kiss_fft_cpx*)fIWPsi.data(), (kiss_fft_cpx*)fIWKPsi.data());

		Scale(fIWKPsi, 1.0 / (1.0*fN*I));

		for (size_t i = 0; i < fN; ++i) {
			Real k = i * Pi / (fDx * fN);
#if 0
			Real s = sin(0.5 * k * fDx);
			Real t = fHbar * (2 * s*s / (fDx*fDx)) / fMass;
#else
			Real t = fHbar * k*k / (2 * fMass);
#endif
			Complex f = exp(-I * (t * fDt));
			fIWKPsi[i] *= f;

			//fKinEnergy += t * fHbar*abs2(fIWKPsi[i]);
			//kinEnergyNorm2 += abs2(fIWKPsi[i]);
		}

		//fKinEnergy /= kinEnergyNorm2;

		fIWKPsi[0] = 0;
		fIWKPsi[fN] = 0;
		for (size_t i = fN + 1; i < 2 * fN; ++i) {
			fIWKPsi[i] = -fIWKPsi[2 * fN - i];
		}

		inv_fft.transform(fIWKPsi.data(), fIWPsi.data());
		//kiss_fft(kiss_cfg_inv, (kiss_fft_cpx*)fIWKPsi.data(), (kiss_fft_cpx*)fIWPsi.data());

		Scale(fIWPsi, 1.0 / (2.0 * I));
		std::copy(fIWPsi.begin(), fIWPsi.begin() + fN, tpsi.begin());

	}
}

void SystemHalfVTHalfV::update_psi()
{
	// exp( (h1+h2£©t) ~ 1 + h1 t + h2 t + 0.5 (h1 + h2 + h1 h2 + h2 h1) t^2
	// exp( h1 t) ~ 1 + h1 t + 0.5 (h1 t)^2
	// exp( h2 t) ~ 1 + h2 t + 0.5 (h2 t)^2
	// exp( 0.5 h1 t) exp( h2 t) exp( 0.5 h1 t) 
	// ~ 1 + h1 t + h2 t + 0.5 (0.5 h1 t)^2 + 0.5 (h2 t)^2 + 0.5 (0.5 h1 t)^2 
	//   + 0.5 h1 h2 t^2 + 0.5 h2 h1 t^2 + ( 0.5  h1 t)^2
	// =  1 + h1 t + h2 t + 0.5 (h1 t)^2 + 0.5 (h2 t)^2 + 0.5 (h1 h2 + h2 h1)t^2 
	// ~ exp( (h1 + h2)t)

	ExpV(fVPsi, fPsi, 0.5);
	ExpT(fTVPsi, fVPsi);
	//Copy(fTVPsi, fVPsi);
	ExpV(fPsi, fTVPsi, 0.5);
	//Copy(fPsi, fVTVPsi);

}

// vpsi = exp(-i/hbar V Dt) psi
void SystemHalfVTHalfV::ExpV(PsiVector &vpsi, PsiVector const &psi, double t)
{
	if (t == 0.5) { // because exp() is very slow
		for (size_t i = 0; i < fN; ++i) {
			vpsi[i] = psi[i] * fExpV0Dot5Dt[i];
		}
	} else {
		Complex f = -1.0 / fHbar * fDt * t;
		for (size_t i = 0; i < fN; ++i) {
			vpsi[i] = psi[i] * exp(f*fV[i] * I);
		}
	}
}

Real SystemHalfVTHalfV::CalKinEn()
{
	if (fBoundaryCondition == BoundaryCondition::Period) {
		//Zero(fFTPsi);

		kissfft<Real> &fft = *(kissfft<Real>*)this->fft_N;
		//kiss_fft_cfg kiss_cfg = kiss_fft_alloc(fN, false, NULL, NULL);

		fft.transform(fPsi.data(), fFTPsi.data());
		//kiss_fft(kiss_cfg, (kiss_fft_cpx*)fPsi.data(), (kiss_fft_cpx*)fFTPsi.data());

		Real kinen = 0;
		Real norm2 = 0;
		for (size_t i = 0; i < fN; ++i) {
			size_t j = i > fN / 2 ? fN - i : i;
			Real p = fHbar * j * 2 * Pi / (fDx * fN);
#if 0
			Real s = sin(0.5 * k * fDx);
			Real t = fHbar * (2 * s*s / (fDx*fDx)) / fMass;
#else
			Real t = p * p / (2 * fMass);
#endif
			kinen += abs2(fFTPsi[i]) * t;
			norm2 += abs2(fFTPsi[i]);
		}
		return kinen / norm2;


	} else if (fBoundaryCondition == BoundaryCondition::InfiniteWall) {

		kissfft<Real> &inv_fft = *(kissfft<Real>*)inv_fft_2N;

		std::copy(fPsi.begin(), fPsi.end(), fIWPsi.begin());
		fIWPsi[0] = 0;
		fIWPsi[fN] = 0;
		for (size_t i = fN + 1; i < 2 * fN; ++i) {
			fIWPsi[i] = -fIWPsi[2 * fN - i];
		}

		inv_fft.transform(fIWPsi.data(), fIWKPsi.data());


		Real kinen = 0;
		Real norm2 = 0;

		for (size_t i = 0; i < fN; ++i) {
			Real p = fHbar * i * Pi / (fDx * fN);
#if 0
			Real s = sin(0.5 * k * fDx);
			Real t = fHbar * (2 * s*s / (fDx*fDx)) / fMass;
#else
			Real t = p * p / (2 * fMass);
#endif
			kinen += abs2(fIWKPsi[i]) * t;
			norm2 += abs2(fIWKPsi[i]);
		}
		return kinen / norm2;

	}
	return 0;
}