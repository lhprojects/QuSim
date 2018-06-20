#define _CRT_SECURE_NO_WARNINGS

#include "System.h"
#include "erfz.h"
#include "Cal.h"

#include "kissfft.hh"

using std::abs;

void dump(PsiVector const &v, char const *fn)
{
	FILE *f = fopen(fn," w");
	for (auto &x : v) {
		fprintf(f, "(% .20lf, % .20lf)\n", x.real(), x.imag());
	}
	fclose(f);
};

void dump(std::vector<Real> const &v, char const *fn)
{
	FILE *f = fopen(fn, " w");
	for (auto &x : v) {
		fprintf(f, "% .20lf\n", x);
	}
	fclose(f);
};

void System::init(char const *psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	char const *vs, Real x0, Real x1, size_t n, BoundaryCondition b,
	Real mass, Real hbar)
{
	testCal();

	fStep = 0;

	fFN = force_normalization;
	fPsiStr = psi;

	fX0 = x0;
	fDx = (x1 - x0) / n;
	fN = n;
	fVStr = vs;
	fBoundaryCondition = b;

	fDt = dt;
	fFNES = force_normalization_each_step;
	fMass = mass;
	fHbar = hbar;

	fPsi.resize(n);
	fVPsi.resize(n);
	fTVPsi.resize(n);
	fVTVPsi.resize(n);

	fV.resize(n);
	initPsi();
	fLastLastPsi = fPsi;
	fLastPsi = fPsi;
	initPotential();

	if (fBoundaryCondition == BoundaryCondition::ExtendInfinity) {
		initFreeParticleProp();
	} else if (fBoundaryCondition == BoundaryCondition::Period) {
		fFTPsi.resize(n);
	} else if (fBoundaryCondition == BoundaryCondition::InfiniteWall) {
		fIWPsi.resize(2 * n);
		fIWKPsi.resize(2 * n);
	}
}

void System::initPsi()
{
	Cal cal(fPsiStr.c_str());
	for (size_t i = 0; i < fN; ++i) {
		Real x = getX(i);
		cal.SetVarVal("x", Complex(x));
		Complex com = cal.Val();
		fPsi[i] = com;
	}

	if (fFN) {
		Scale(fPsi, 1 / sqrt(Norm2()));
	}
	dump(fPsi, "psi.txt");
}

void System::initPotential()
{
	Cal cal(fVStr.c_str());
	for (size_t i = 0; i < fN; ++i) {
		Real x = getX(i);
		cal.SetVarVal("x", Complex(x));
		Complex com = cal.Val();
		fV[i] = com.real();
	}
	//dump(fV, "V.txt");
}

void System::initFreeParticleProp()
{
	fProp.resize(fN);
	for (size_t i = 0; i < fN; ++i) {
		Complex T = fHbar / fMass * fDt;
		double x = i * fDx;
		double K = 2 * Pi / fDx;
		Complex a1 = 0.25*(1.0 + I)*(K*T - 2 * x) / sqrt(T);
		Complex a2 = 0.25*(1.0 + I)*(K*T + 2 * x) / sqrt(T);
		Complex corr = 0.5 *(erfz(a1) + erfz(a2));
		fProp[i] = 0.5*(1.0 - I)*exp(0.5*x*x / T * I) / sqrt(Pi*T)*corr*fDx;
	}
	double sum = 0;
	for (size_t i = 0; i < fN; ++i) {
		sum += (i != 0 ? 2 : 1) * abs(fProp[i])*abs(fProp[i]);
	}
	//printf("prop %.20lf\n", sum);
}



void System::ExpT(PsiVector &tpsi, PsiVector const &psi)
{
	Zero(tpsi);
	if (fBoundaryCondition == BoundaryCondition::ExtendInfinity) {
		for (size_t i = 0; i < fN; ++i) {
			for (size_t j = 0; j < fN; ++j) {
				tpsi[i] += psi[j] * fProp[(size_t)abs((Int)i - (Int)j)];
			}
		}
	} else if (fBoundaryCondition == BoundaryCondition::Period) {
		//Zero(fFTPsi);

		kissfft<Real> fft(fN, false);
		kissfft<Real> inv_fft(fN, true);
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
			Real t = fHbar * (2*s*s/(fDx*fDx)) / fMass;
#else
			Real t = fHbar * k*k / (2 * fMass);
#endif
			Complex f = exp(-I * (t * fDt));
			fFTPsi[i] *= f;
		}

		inv_fft.transform(fFTPsi.data(), tpsi.data());
		//kiss_fft(kiss_cfg_inv, (kiss_fft_cpx*)fFTPsi.data(), (kiss_fft_cpx*)tpsi.data());

		Scale(tpsi, 1.0 / sqrt(1.0*fN));

		//free(kiss_cfg);
		//free(kiss_cfg_inv);


	} else if (fBoundaryCondition == BoundaryCondition::InfiniteWall) {
	
		kissfft<Real> inv_fft(2*fN, true);
		//kiss_fft_cfg kiss_cfg_inv = kiss_fft_alloc(2*fN, true, NULL, NULL);

		std::copy(psi.begin(), psi.end(), fIWPsi.begin());
		fIWPsi[0] = 0;
		fIWPsi[fN] = 0;
		for (size_t i = fN + 1; i < 2 * fN; ++i) {
			fIWPsi[i] = -fIWPsi[2*fN - i];
		}

		inv_fft.transform(fIWPsi.data(), fIWKPsi.data());
		//kiss_fft(kiss_cfg_inv, (kiss_fft_cpx*)fIWPsi.data(), (kiss_fft_cpx*)fIWKPsi.data());

		Scale(fIWKPsi, 1.0/(1.0*fN*I));

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
		}
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

void System::step()
{
	// exp( (h1+h2£©t) ~ 1 + h1 t + h2 t + 0.5 (h1 + h2 + h1 h2 + h2 h1) t^2
	// exp( h1 t) ~ 1 + h1 t + 0.5 (h1 t)^2
	// exp( h2 t) ~ 1 + h2 t + 0.5 (h2 t)^2
	// exp( 0.5 h1 t) exp( h2 t) exp( 0.5 h1 t) 
	// ~ 1 + h1 t + h2 t + 0.5 (0.5 h1 t)^2 + 0.5 (h2 t)^2 + 0.5 (0.5 h1 t)^2 
	//   + 0.5 h1 h2 t^2 + 0.5 h2 h1 t^2 + ( 0.5  h1 t)^2
	// =  1 + h1 t + h2 t + 0.5 (h1 t)^2 + 0.5 (h2 t)^2 + 0.5 (h1 h2 + h2 h1)t^2 
	// ~ exp( (h1 + h2)t)

	fLastLastPsi = fLastPsi;
	fLastPsi = fPsi;
	ExpV(fVPsi, fPsi, 0.5);
	ExpT(fTVPsi, fVPsi);
	//Copy(fTVPsi, fVPsi);
	ExpV(fVTVPsi, fTVPsi, 0.5);

	Copy(fPsi, fVTVPsi);
	if (fFNES) {
		Scale(fPsi, 1.0 / sqrt(Norm2()));
	}

	++fStep;
}

Real System::KinEn()
{
	if (fBoundaryCondition == BoundaryCondition::ExtendInfinity) {
		return 0;
	} else if (fBoundaryCondition == BoundaryCondition::Period) {
		//Zero(fFTPsi);

		kissfft<Real> fft(fN, false);
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
			Real t = p*p / (2 * fMass);
#endif
			kinen += abs2(fFTPsi[i]) * t;
			norm2 += abs2(fFTPsi[i]);
		}
		return kinen / norm2;


	} else if (fBoundaryCondition == BoundaryCondition::InfiniteWall) {

		kissfft<Real> inv_fft(2 * fN, true);

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
			Real t = p*p / (2 * fMass);
#endif
			kinen += abs2(fIWKPsi[i]) * t;
			norm2 += abs2(fIWKPsi[i]);
		}
		return kinen / norm2;

	}
	return 0;
}

Real System::EnPartialT()
{
	Complex en = 0;
	for (size_t i = 0; i < fN; ++i) {
		en += I*fHbar*conj(fPsi[i])*(fPsi[i] - fLastLastPsi[i])/(2.0*fDt)*fDx;
	}
	return en.real();
}
