#define _CRT_SECURE_NO_WARNINGS

#include "System.h"
#include "erfz.h"
#include "Cal.h"
#include "kissfft.hh"

using std::abs;

void System::init(char const *psi, bool fn, Real dt,
	char const *vs, Real x0, Real x1, size_t n, BoundaryCondition b,
	Real mass, Real hbar)
{
	testCal();

	fStep = 0;

	fFN = fn;
	fPsiStr = psi;

	fX0 = x0;
	fDx = (x1 - x0) / (n - 1);
	fN = n;
	fVStr = vs;
	fBoundaryCondition = b;

	fDt = dt;
	fMass = mass;
	fHbar = hbar;

	fPsi.resize(n);
	fVPsi.resize(n);
	fTVPsi.resize(n);
	fVTVPsi.resize(n);

	fV.resize(n);
	initPsi();
	initPotential();

	if (fBoundaryCondition == BoundaryCondition::ExtendInfinity) {
		initFreeParticleProp();
	} else if (fBoundaryCondition == BoundaryCondition::Period) {
		fFTPsi.resize(n);
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
		Real n = 1 / sqrt(Norm2());
		for (size_t i = 0; i < fN; ++i) {
			fPsi[i] *= n;
		}
	}
}

void System::initPotential()
{
	Cal cal(fVStr.c_str());
	for (size_t i = 0; i < fN; ++i) {
		Real x = getX(i);
		cal.SetVarVal("x", Complex(x));
		Complex com = cal.Val();
		fV[i] = abs(com);
	}
}

void System::initFreeParticleProp()
{
	fProp.resize(fN);
	for (size_t i = 0; i < fN; ++i) {
		double T = fHbar / fMass * fDt;
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

		fft.transform(fPsi.data(), fFTPsi.data());
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
		Scale(tpsi, 1.0 / sqrt(1.0*fN));

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

	

	ExpV(fVPsi, fPsi, 0.5);
	ExpT(fTVPsi, fVPsi);
	//Copy(fTVPsi, fVPsi);
	ExpV(fVTVPsi, fTVPsi, 0.5);

	Copy(fPsi, fVTVPsi);

	++fStep;
}
