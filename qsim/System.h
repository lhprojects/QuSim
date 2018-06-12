#pragma once



#include "erfz.h"
#include "Cal.h"
#include <vector>
#include <complex>
#include <cmath>
#include <functional>


using Real = double;
using Complex = std::complex<Real>;
using PsiVector = std::vector<Complex>;
using Int = int64_t;
using UInt = uint64_t;
static constexpr Complex I = Complex(0, 1);
static constexpr Real Pi = 3.141592653589793238462643383279502884197169399375;

#include <cmath>

inline Real abs2(Complex const &c)
{
	return real(c)*real(c) + imag(c)*imag(c);
}

struct System {


	Int fStep;
	Real fDt;
	Real fMass;
	Real fHbar;
	Real fDx;
	Real fX0;
	size_t fN;
	std::string fVStr;
	std::string fPsiStr;
	bool fFN;
	std::vector<Complex> fPsi;
	std::vector<Complex> fVPsi;
	std::vector<Complex> fTVPsi;
	std::vector<Complex> fVTVPsi;

	std::vector<Complex> fProp;
	std::vector<Real> fV;

	System()
	{
		fN = 0;
	}
	Real getX(size_t i)
	{
		return fDx * i + fX0;
	}

	void init(char const *psis, bool fn, Real dt, Real x0, Real x1, char const *vs, size_t n, Real mass, Real hbar = 1)
	{
		testCal();

		fStep = 0;
		fPsiStr = psis;
		fVStr = vs;
		fFN = fn;
		fDt = dt;
		fDx = (x1 - x0) / n;
		fX0 = x0;
		fN = n;
		fMass = mass;
		fHbar = hbar;
		fPsi.resize(n);
		fVPsi.resize(n);
		fTVPsi.resize(n);
		fVTVPsi.resize(n);

		fProp.resize(n);
		fV.resize(n);
		initPsi();
		initPotential();
		initFreeParticleProp();
	}


	void initPsi()
	{
		Cal cal(fPsiStr.c_str());
		for (size_t i = 0; i < fN; ++i) {
			Real x = getX(i);
			cal.SetVarVal("x", Complex(x));
			Complex com = cal.Val();
			fPsi[i] = com;
		}
		if (fFN) {
			Real n = 1 / Norm();
			for (size_t i = 0; i < fN; ++i) {
				fPsi[i] *= n;
			}
		}
	}

	void initPotential()
	{
		Cal cal(fVStr.c_str());
		for (size_t i = 0; i < fN; ++i) {
			Real x = getX(i);
			cal.SetVarVal("x", Complex(x));
			Complex com = cal.Val();
			fV[i] = abs(com);
		}
	}

	void initFreeParticleProp()
	{
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

	// vpsi = exp(-i/hbar V Dt) psi
	void ExpV(PsiVector &vpsi, PsiVector const &psi, double t)
	{
		Real f = -1.0 / fHbar * fDt * t;
		for (size_t i = 0; i < fN; ++i) {
			vpsi[i] = psi[i] * exp(f*fV[i] * I);
		}
	}

	void ExpT(PsiVector &tpsi, PsiVector const &psi)
	{
		Zero(tpsi);
		for (size_t i = 0; i < fN; ++i) {
			for (size_t j = 0; j < fN; ++j) {
				tpsi[i] += psi[j] * fProp[(size_t)abs((Int)i - (Int)j)];
			}
		}
	}

	void Zero(PsiVector &psi)
	{
		for (size_t i = 0; i < fN; ++i) {
			psi[i] = 0.;
		}
	}

	void Add(PsiVector &psi, PsiVector const &psi1, PsiVector const &psi2)
	{
		for (size_t i = 0; i < fN; ++i) {
			psi[i] = psi1[i] + psi2[i];
		}
	}

	void Copy(PsiVector &psi, PsiVector const &psi1)
	{
		for (size_t i = 0; i < fN; ++i) {
			psi[i] = psi1[i];
		}
	}

	void step()
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

	Real Time()
	{
		return fStep * fDt;
	}
	Real Xavg()
	{
		Real norm2 = 0;
		for (size_t i = 0; i < fN; ++i) {
			norm2 += abs2(fPsi[i])*getX(i)*fDx;
		}
		return norm2/Norm2();
	}

	Real PotEn()
	{
		Real norm2 = 0;
		for (size_t i = 0; i < fN; ++i) {
			norm2 += abs(fPsi[i])*abs(fPsi[i])*fV[i];
		}
		return norm2 / Norm2();
	}

	Real KinEn()
	{
		Real norm2 = 0;
		for (size_t i = 1; i < fN - 1; ++i) {
			norm2 += abs((fPsi[i - 1] / fPsi[i] + fPsi[i + 1] / fPsi[i] - 2.0)*fHbar / (2 * fMass));
		}
		return norm2 / Norm2();
	}

	Real Norm2()
	{
		Real norm2 = 0;
		for (size_t i = 0; i < fN; ++i) {
			norm2 += abs2(fPsi[i])*fDx;
		}
		return norm2;
	}

	Real Norm()
	{
		return sqrt(Norm2());
	}

};

